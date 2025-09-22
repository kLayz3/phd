#pragma once

#include "libs.hh"
#include "TTree.h"
#include "TChain.h"
#include <atomic>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <optional>
#include "TAnalysisWorker.hxx"
#include "TContainer.h"
#include "dbg.hh"

/* Passing in this define to purposefully undef
 * the singlethreaded-ness. */
#if defined(ANALYSIS_MULTITHREADED)
#	undef ANALYSIS_SINGLETHREADED
#endif

/* Passed from the build tool. */
#if !defined(ANALYSIS_SINGLETHREADED)
#	if !defined(POOL_MAX_THREADS_)
#		define POOL_MAX_THREADS_ 10
#	endif
#else
#	warning "Running single-threaded. Possibly slower for complex `ProcessEntry` calls."
#	define POOL_MAX_THREADS 100 /* Dummy value. */
#endif

namespace util {
	template<typename T, typename = void> struct worker_of {};
	template<typename T> struct worker_of<TAnalysisWorker<T>, void> { using type = T; };

	template<typename Tuple, typename F, std::size_t... Is>
	void _for_each_in_tuple_impl(Tuple&& t, F&& f, std::index_sequence<Is...>) {
		(f (std::get<Is>(std::forward<Tuple>(t))), ...);
	}

	template<typename Tuple, typename F>
	void for_each_in_tuple(Tuple&& t, F&& f) {
		constexpr std::size_t N = std::tuple_size_v<std::decay_t<Tuple>>;
		_for_each_in_tuple_impl(std::forward<Tuple>(t), std::forward<F>(f), 
			std::make_index_sequence<N>{});
	}
}

template<typename... Ts>
class TAnalysisPool final {
	static_assert(sizeof...(Ts) <= POOL_MAX_THREADS_, 
		"TAnalysisPool template instantiated with over-the-top capacity. Consider lowering.");
	static_assert((std::is_base_of_v<TProcessor, Ts> && ...), "All the (tuple) base workers must inherit from TProcessor!");
public:
	static constexpr size_t Size() { return sizeof...(Ts); }
	
private:
	std::tuple<
		TAnalysisWorker<Ts>...
	> _pool;

	bool _stop = true;     // Flagging this will join the workers' threads back to the main
	bool _is_valid = true; // Validity flag. Moved objects will be marked as invalid.
	
	template<size_t I>
	using worker_t = std::tuple_element_t<I, decltype(_pool)>;
	template<size_t I>
	using base_t = typename util::worker_of<worker_t<I>>::type;

	template<size_t I>
	static TProcessor* get_at(TAnalysisPool* self) noexcept {
		return static_cast<TProcessor*>(&std::get<I>(self->_pool));
	}
	
	using Getter = TProcessor* (*)(TAnalysisPool*);

	template<std::size_t... Is>
	static constexpr std::array<Getter, sizeof...(Is)>
	make_getter_table(std::index_sequence<Is...>) noexcept { return { &get_at<Is>... }; }

public:
	// TTree objects are owned and handled by gROOT. Here we expose just the raw handles.
	std::variant<TTree*, TChain*> in;
	std::optional<TTree*> out;

	TAnalysisPool() = default;
	explicit TAnalysisPool(std::tuple<TAnalysisWorker<Ts>...> w) : _pool(std::move(w)) {}

	TAnalysisPool(const TAnalysisPool& )            = delete;
	TAnalysisPool& operator=(const TAnalysisPool& ) = delete;
	TAnalysisPool(TAnalysisPool&& ) noexcept            = default;
	TAnalysisPool& operator=(TAnalysisPool&& ) noexcept = default;

	~TAnalysisPool() = default;

	/* Hacky part. Adding a new worker changes the type of the object. */

	// Reference-qualified: take an l/rvalue ref of an existing worker and
	// eat up the object. Now we own it. 
	// Moves *this* pool into a new one.
	template<typename W,
		typename Decayed = std::decay_t<W>,
		typename U = typename util::worker_of<Decayed>::type
	> auto push_worker(W&& w) && -> TAnalysisPool<Ts..., U> {
		auto new_pool = std::tuple_cat(
			std::move(_pool),
			std::tuple<TAnalysisWorker<U>>(std::move(w))
		);
		this->_is_valid = false; // Invalidate the current object.
		return TAnalysisPool<Ts..., U> ( std::move(new_pool) );
	}
	
	// Emplace-style construction, again move *this* pool into a new one.
	template<typename U, typename... Args,
		typename std::enable_if<
			std::is_constructible_v<
				TAnalysisWorker<U>,
				Args...
			>
		>::type* = nullptr
	> auto emplace_worker(Args&&... args) && -> TAnalysisPool<Ts..., U> {
		auto new_pool = std::tuple_cat(
			std::move(_pool),
			std::make_tuple(TAnalysisWorker<U>(std::forward<Args>(args)... ))
		);
		return TAnalysisPool<Ts..., U>(std::move(new_pool));
	}

	bool IsStopped() const noexcept { return _stop; }

	/**
	 * Returns the pointer 'T*' of the object wrapped up in the 'TAnalysisWorker<T>',
	 * at the tuple index 'I'
	 */
	template<std::size_t I>
	auto* GetWorker() noexcept {
		static_assert(I < Size(), "Request for worker index outside of the tuple size.");
		return static_cast<base_t<I>*>(&std::get<I>(_pool));
	}
	template<std::size_t I>
	auto* GetWorker() const noexcept {
		static_assert(I < Size(), "Request for worker index outside of the tuple size.");
		return static_cast<const base_t<I>*>(&std::get<I>(_pool));
	}
	// Runtime version, can throw.
	TProcessor* GetWorker(size_t i) {
		if(i >= Size()) ERROR("Request for worker index outside of the tuple size.");
		static constexpr auto table = make_getter_table(std::make_index_sequence<Size()>{});

		return table[i](this);
	}

	Int_t GetEntry(Long64_t entry, Int_t getall = 0) const { 
		if(in.valueless_by_exception()) 
			ERROR("Valueless input TTree/TChain. Invalid");
		else return std::visit([=](const auto& obj) { return obj->GetEntry(entry,getall); }, this->in); 
	}

	Long64_t GetEntries() const noexcept { 
		if(std::holds_alternative<TTree*>(in)) 
			return std::get<TTree*>(in)->GetEntries();
		else if(std::holds_alternative<TChain*>(in)) 
			return std::get<TChain*>(in)->GetEntries();
		else return 0;
	}

	/**
	 * Calls the internal `Fill` of the (optional) TTree* output object
	 */
	Int_t FillOutput() const noexcept { 
		if(out and *out) return (*out)->Fill();
		else return 0;
	}
	
	/**
	 * Start the underlying thread of every worker.
	 */
	void Start() {
		if(!_is_valid) ERROR("Called in invalidated object.");
		// State prior to this call must be stopped.
		if(!_stop) ERROR("Must be in a stopped state before calling start.");
		_stop = false;
#if !defined(ANALYSIS_SINGLETHREADED)
		util::for_each_in_tuple(_pool, [](auto& w) {
			w._stop.store(false, std::memory_order_release);
			w.Start();
		});
#endif
	}
	
	/**
	 * Sequentially calls `T::ProcessEntry()` of each underlying worker. 
	 * */
	void ProcessEntry() noexcept {
		std::apply([](auto&... ws) { (..., ws.ProcessEntry()); }, _pool);
	}

	/**
	 * Mark the flag to wake-up for all workers simultaneously (multithreaded mode).
	 * Call underlying workers' ProcessEntry() sequentially in single threaded mode. 
	 */
	void AssignWork() noexcept {
#if defined(ANALYSIS_SINGLETHREADED)
		this->ProcessEntry();
#else
		std::apply([](auto&... ws) { 
				(..., ws._has_work.store(true, std::memory_order_release)); 
			}, _pool); 
#endif
	};

	/**
	 * Block the main thread until the execution of all the workers is marked finished. 
	 * No-op for single threaded mode. 
	 */
	void Await() const noexcept {
#if !defined(ANALYSIS_SINGLETHREADED)
		while(true) {
			bool all_done = std::apply([](auto&... ws) {
					return (...&& ws.IsDone()); 	
				}, _pool);
			if(all_done) break;
#if defined(__x86_64__)
			_mm_pause();
#endif
		}
#endif
	};

	/**
	 * Stops execution of all the workers and joins their internal threads.
	 * No-op for single threaded mode. 
	 */
	void Stop() {
		if(!_is_valid) ERROR("Called in invalidated object.");
		_stop = true;
		
		util::for_each_in_tuple(_pool, [](auto& w) {
			if(w._thread.joinable()) {
				w._stop.store(true, std::memory_order_release);
				w._thread.join();
			}
		});	
	};

	/**
	 * At the end, write the objects and the output TTree into the TFile.
	 * Note that this will make the pool instance unusable for the remainder of the program.
	 * Returns number of bytes written across all the output containers.
	 */
	Int_t Write() {
		if(!_is_valid) ERROR("Called in invalidated object.");
		Int_t r = std::apply([](auto&... ws) {
				return (...+ ws.Write());
			}, _pool
		);
		// Tree handle just gets written into its directory.
		// Check for possible mishandling though.
		if(out) {
			if(!(*out)) ERROR("Output TTree handle given, but is nullptr? Cannot write.");

			TDirectory* dir = (*out)->GetDirectory();
			if(!dir) ERROR("Output TTree handle given, but is not attached to any directory.");

			TFile *f = dynamic_cast<TFile*>(dir);
			if(!f) ERROR("Output TTree handle belongs to a non-file directory.");
			if(!f->IsWritable()) ERROR("Output TTree's associated file is not writable. Inside 'main()' define it after opening an output file, or call its 'SetDirectory()` method.");

			r += (*out)->Write();
			dbg("Successfully written the output TTree.");
		} else {
			WARN("Calling Write but output TTree* not specified. Is OK.");
		}
		_is_valid = false;
		return r;
	};
	
}; // TAnalysisPool
