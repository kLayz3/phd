#pragma once

#include "libs.hh"
#include "TTree.h"
#include "TChain.h"
#include <atomic>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include "AuxFunctions.hh"
#include <variant>
#include "TAnalysisWorker.hxx"
#include "TContainer.hxx"
#include "TFile.h"
#include "TKey.h"
#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"

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

	// Pool indirectly owns all the TFile handles.
	std::string in_file{}; // can be empty.
	std::string out_file{}; // can be empty.
	std::string out_rntuple{}; // can be empty

public:
	using Empty = std::monostate;
	std::variant<
		Empty,
		TTree*,  // non-owning
		TChain*, // non-owning
		std::unique_ptr<ROOT::Experimental::RNTupleReader>
	> _in_reader;
	std::unique_ptr<ROOT::Experimental::RNTupleWriter> _out_writer; // can be nullptr.

	TAnalysisPool() = default;
	TAnalysisPool(
		std::variant<std::string, TTree*, TChain*> input_handle, 
		std::string out_fname,
		std::string out_rntuple_name
	) : TAnalysisPool() {
		SetInput(input_handle), SetOutput(out_fname, out_rntuple_name); 
		// No need to perfectly forward here. This call gets executed only once pretty much. 
	}
	explicit TAnalysisPool(
		std::tuple<TAnalysisWorker<Ts>...>&& w,
		std::string&& in_file,
		std::string&& out_file,
		std::string&& out_rntuple,
		std::variant<Empty, TTree*, TChain*,
			std::unique_ptr<ROOT::Experimental::RNTupleReader>
		>&& _in_reader,
		std::unique_ptr<ROOT::Experimental::RNTupleWriter>&& _out_writer
	) : _pool(std::move(w)), in_file(std::move(in_file)), 
		out_file(std::move(out_file)), out_rntuple(std::move(out_rntuple)),
		_in_reader(std::move(_in_reader)), _out_writer(std::move(_out_writer)) 
	{}

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
		return TAnalysisPool<Ts..., U> ( 
			std::move(new_pool),
			std::move(in_file), std::move(out_file), std::move(out_rntuple),
			std::move(_in_reader), std::move(_out_writer)
		);
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
		return TAnalysisPool<Ts..., U>(
			std::move(new_pool),
			std::move(in_file), std::move(out_file), std::move(out_rntuple),
			std::move(_in_reader), std::move(_out_writer)
		);
	}
	
	void SetInput(std::variant<std::string, TTree*, TChain*> input) {
		if(input.valueless_by_exception())
			ERROR("Passed valueless variant?\n");
		else if(std::holds_alternative<TTree*>(input))
			_in_reader = std::get<TTree*>(input);
		else if(std::holds_alternative<TChain*>(input))
			_in_reader = std::get<TChain*>(input);
		else {
			std::string file_name = std::get<std::string>(input);
			if(!util::is_file_readable(file_name))
				ERROR("Unable to read a ROOT file: " EMPH(%s\n), file_name.c_str());
			
			/* RNTuple must exist in this file. 
			 * Find the first one and mark it. Multiple RNTuple objects qualify
			 * undefined behaviour. */
			auto f = std::make_unique<TFile>(file_name.c_str(), "READ");	
			if(f->IsZombie() || !f->IsOpen())
				ERROR("Able to open a ROOT file, but unable to make a TFile hook in: " EMPH(%s\n), file_name.c_str());
			
			std::string rntuple_name;

			for(TObject* _k : *f->GetListOfKeys()) {
				TKey* k = dynamic_cast<TKey*>(_k);
				if(!k) continue;
				TClass* cl = TClass::GetClass(k->GetClassName());
				if(cl && cl->InheritsFrom(ROOT::RNTuple::Class())) {
					rntuple_name = std::string(k->GetName());
					break;
				}
			}
			if(rntuple_name.empty())
				ERROR("Able to cleanly read the ROOT file, but couldn't find the ROOT::RNTuple object in: " EMPH(%s\n), file_name.c_str());
			{ std::unique_ptr<TFile> _ = std::move(f); } // Release the resource.

			this->in_file = file_name;
			_in_reader = ROOT::Experimental::RNTupleReader::Open(std::move(TContainerBase::ReleaseModelRead()), rntuple_name, file_name);
		}
	}

	void SetOutput(std::string file_name, std::string rntuple_name) {
		if(file_name.empty() || rntuple_name.empty())
			ERROR("Argument strings empty? \'%s\' , \'%s\'\n", file_name.c_str(), rntuple_name.c_str());

		this->out_file = std::move(file_name);
		this->out_rntuple = std::move(rntuple_name);

		_out_writer = ROOT::Experimental::RNTupleWriter::Recreate(std::move(TContainerBase::ReleaseModelWrite()), out_rntuple, out_file);
		if(!_out_writer)
			ERROR("Output file unable to be (re)created. Is the file path valid: " EMPH(%s\n), out_file.c_str());
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
		if(_in_reader.valueless_by_exception()) 
			ERROR("Valueless input TTree/TChain/RNTuple. Invalid");
		else if(std::holds_alternative<Empty>(_in_reader)) 
			ERROR("Empty input TTree/TChain/RNTuple. Invalid");
		else if(std::holds_alternative<TTree*>(_in_reader)) 
			return std::get<TTree*>(_in_reader)->GetEntry(entry, getall);
		else if(std::holds_alternative<TChain*>(_in_reader)) 
			return std::get<TChain*>(_in_reader)->GetEntry(entry, getall);
		else {
			std::get<
				std::unique_ptr<ROOT::Experimental::RNTupleReader>
			> (_in_reader)->LoadEntry(entry);
			return 0;
		}
	}

	Long64_t GetEntries() const { 
		if(_in_reader.valueless_by_exception()) 
			ERROR("Valueless input TTree/TChain/RNTuple. Invalid");
		else if(std::holds_alternative<Empty>(_in_reader)) 
			ERROR("Empty input TTree/TChain/RNTuple. Invalid");
		else if(std::holds_alternative<TTree*>(_in_reader)) 
			return std::get<TTree*>(_in_reader)->GetEntries();
		else if(std::holds_alternative<TChain*>(_in_reader)) 
			return std::get<TChain*>(_in_reader)->GetEntries();
		else 
			return static_cast<Long64_t>( std::get<
				std::unique_ptr<ROOT::Experimental::RNTupleReader>
			> (_in_reader)->GetNEntries());
	}

	/**
	 * Calls the internal `Fill` of the (optional) RDTuple output object
	 */
	Int_t Fill() const noexcept { 
		if(_out_writer) return _out_writer->Fill();
		else return 0;
	}
	
	/**
	 * Start the underlying thread of every worker.
	 */
	void Start() {
		if(!_is_valid) 
			ERROR("Called in invalidated object.");
		// State prior to this call must be stopped.
		if(!_stop) 
			ERROR("Must be in a stopped state before calling start.");
		if(_in_reader.valueless_by_exception()) 
			ERROR("Valueless input TTree/TChain/RNTuple. Invalid");
	
		/* Check if the reader state isn't Empty, and the pointer variant valid. */
		if(std::visit([](const auto& v) { 
			using X = std::decay_t<decltype(v)>;
			if constexpr(std::is_same_v<X, Empty>)
				return true;
			else
				return !v; 
		}, _in_reader ) ) {
			ERROR("Reader input is in valued, but either empty or null state?\n");
		}

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
	 * At the end, write the objects and the output RNTuple into the TFile.
	 * Note that this will make the pool instance unusable for the remainder of the program.
	 * Returns number of bytes written across all the output containers.
	 * Since the RNTuple object owns the file, it needs to write first.
	 */
	Int_t Write() {
		if(out_file.empty() or !_is_valid) {
			WARN("Called Write but pool object is either invalid or output string not specified.\n");
			return 0;
		}
	
		/* First write the RNTuple. */
		_out_writer->CommitCluster();
		_out_writer.reset();

		/* Tail calls of following sequence:
		 * `TAnalysisWorker<Whatever>::Write` calls
		 * `Whatever::Write` calls `TContainer<???>::Write`.
		 * Now, this **should** write to the same file as RNTupleWriter, or it could be
		 * a different file, but then weird things happen. */

		Int_t r = std::apply([](auto&... ws) {
				return (...+ ws.Write());
			}, _pool
		);

		WARN("Successfully written %d bytes the output file.\n", r);
		_is_valid = false;
		return r;
	};
	
}; // TAnalysisPool
