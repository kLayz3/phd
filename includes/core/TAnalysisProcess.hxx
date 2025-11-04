#pragma once

#include "TContainer.hxx"
#include "libs.hh"
#include "TTree.h"
#include "TChain.h"
#include <atomic>
#include <boost/lockfree/policies.hpp>
#include <cstring>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include "AuxFunctions.hh"
#include <variant>
#include "TProcessor.hxx"
#include "TFile.h"
#include "TKey.h"
#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <ROOT/RNTupleReader.hxx>

#include <atomic>
#include <thread>
#include "boost/lockfree/spsc_queue.hpp"

#if __has_include(<immintrin.h>)
#   define __HAS_SMALL_INTEL_SPIN
#	include <immintrin.h>
#endif

#define N_EVENTS_PER_BATCH 512;

/* For ROOT 6.36+, change this to ROOT instead of `ROOT::Experimental`. */
template<u32, typename...> class TAnalysisPool;
namespace RExp = ROOT::Experimental;

namespace util {

template<typename T, typename = void>
struct has_process_entry : std::false_type {};

template<typename T>
struct has_process_entry<T, std::void_t<decltype(std::declval<T&>().ProcessEntry())>> 
	: std::is_same<
		decltype(std::declval<T&>().ProcessEntry()),
		void
	> {};

struct Job { u64 first, last; };
using JobQueue = boost::lockfree::spsc_queue <
	Job, boost::lockfree::capacity<128>
>;

struct IOInfo {
	struct {
		std::string fname;
	} in;
	struct {
		std::string fname;
		std::string out_rnname;
	} out;
};

struct RNPerThreadReader {
	std::unique_ptr<RExp::RNTupleModel>  _model;
	std::unique_ptr<RExp::RNTupleReader> _reader;
};

struct TTreePerThreadReader {
	std::unique_ptr<TFile> _file;
	TTree* _tree;
};

using PerThreadReader = Variant <
	RNPerThreadReader,
	TTreePerThreadReader
>;

/** 
 * Return codes:
 * 0  => properly set-up.
 * 1  => variant in empty state
 * 2  => switched to TTree variant but file handle is null. 
 * 4  => switched to TTree variant but ttree handle is null. 
 * 8  => switched to RNTuple variant but model isn't null.
 * 16 => switched to RNTuple variant but reader is null.
 */
inline int GetValidity(const PerThreadReader& reader) {
	int r = 0;
	if(IsEmpty(reader))   r |= 0x1;
	else if(auto* p = std::get_if<TTreePerThreadReader>( &reader);
		p != nullptr) {
		if(! p->_file )   r |= 0x2;
		if(! p->_tree )   r |= 0x4;
	}
	else if(auto* p = std::get_if<RNPerThreadReader>( &reader);
		p != nullptr) {
		if( p->_model )   r |= 0x8;
		if(! p->_reader ) r |= 0x10;
	}
	return r;
}

struct PerThreadWriter {
	/* Thread [0]'s instance owns the parallelwriter, other threads
	 * will default to be holding only the raw pointer to it. */
	Variant <
		std::unique_ptr<RExp::RNTupleParallelWriter>,
		RExp::RNTupleParallelWriter*
	> pwriter;
	std::shared_ptr<RExp::RNTupleFillContext> ctx;
	std::unique_ptr<RExp::REntry> entry;

	void Reset() {
		entry.reset();
		ctx.reset();
		auto* p = std::get_if<std::unique_ptr<RExp::RNTupleParallelWriter>>(&pwriter);
		if(p != nullptr)
			p->reset();
	}

	/** 
	 * Return codes:
	 * 0 => properly set-up.
	 * 1 => pwriter variant in empty state
	 * 2 => pwriter non-empty state, but is null.
	 * 4 => context in null state.
	 * 8 => entry in null state.
	 */
	int GetValidity() const {
		int r = 0;
		if(IsEmpty(pwriter)) r |= 0x1;
		else if ( /* Try to cast to 1st type instance. */
			auto* p = std::get_if<std::unique_ptr<RExp::RNTupleParallelWriter>>(&pwriter);
			p != nullptr && *p == nullptr 
		)                          r |= 0x2;
		else if ( /* Try to cast to 2nd type instance. */
			auto* p = std::get_if<RExp::RNTupleParallelWriter*>(&pwriter);
			p != nullptr && *p == nullptr 
		)                          r |= 0x2;
		if(! ctx)                  r |= 0x4;
		if(! entry)                r |= 0x8;
		return r;
	}
};
} // namespace util

/**
 * Represents the full analysis process, where an input entry 
 * will go through all of the `T::ProcessEntry()` for each of the 
 * sequential processes to create a complete output event.
 *
 * It holds a queue of entries that the main thread distributes, and encapsulates
 * all the resources a worker thread will need.
 */
template<typename... Ts>
struct alignas(util::CL) TAnalysisProcess final {
	static_assert((util::is_base_of_template<TProcessor, Ts>::value && ...),
		"All the inderlying subprocess types <Ts> must inherit from TProcessor<...>.");
	static_assert((util::has_process_entry<Ts>::value && ...),
		"All the inderlying subprocess types <Ts> need a public `void ProcessEntry() noexcept` method implemented.");
	static_assert((std::is_copy_constructible_v<Ts> && ...),
		"All the inderlying subprocess types <Ts> need a copy ctor.");
	static_assert((std::is_copy_assignable_v<Ts> && ...),
		"All the inderlying subprocess types <Ts> need copy assignment op.");
	static_assert((std::is_move_constructible_v<Ts> && ...),
		"All the inderlying subprocess types <Ts> need a move ctor.");
	static_assert((std::is_move_assignable_v<Ts> && ...),
		"All the inderlying subprocess types <Ts> need move assignment op.");

	template<u32, typename...> friend class TAnalysisPool;
	static constexpr size_t Size() { return sizeof...(Ts); }

private:
	util::JobQueue q;
	std::atomic<bool> _running {false}; // Flagging this will join the workers' threads back to the main
	bool _do_write {true};
	std::thread _thread;

	std::tuple<Ts...> _proc;

	util::IOInfo info;
	util::PerThreadReader reader;
	util::PerThreadWriter writer;
	
private:
	template<size_t I>
	using base_t = std::tuple_element_t<I, decltype(_proc)>;

	template<size_t I>
	static TProcessorBase* get_at(TAnalysisProcess* self) noexcept {
		return static_cast<TProcessorBase*>(&std::get<I>(self->_proc));
	}
	
	using Getter = TProcessorBase* (*)(TAnalysisProcess*);

	template<std::size_t... Is>
	static constexpr std::array<Getter, sizeof...(Is)>
	make_getter_table(std::index_sequence<Is...>) noexcept { return { &get_at<Is>... }; }

public:
	TAnalysisProcess() = default;
	TAnalysisProcess(std::string file_in, std::string file_out, std::string rn_out) :
		TAnalysisProcess() {
			info.in.fname  = std::move(file_in);
			info.out.fname = std::move(file_out);
			info.out.out_rnname = std::move(rn_out);
		}
	
	explicit TAnalysisProcess(
		std::tuple<Ts...>&& w,
		util::IOInfo&& _info,
		util::PerThreadReader&& _reader,
		util::PerThreadWriter&& _writer
	) : _proc(std::move(w)), info(std::move(_info)), 
		reader(std::move(_reader)), writer(std::move(_writer)) {}

	TAnalysisProcess(const TAnalysisProcess& rhs) : TAnalysisProcess{} {
		rhs.Clone(*this);
	}

	TAnalysisProcess& operator=(const TAnalysisProcess& rhs) {
		rhs.Clone(*this);
		return *this;
	}

	/* Move ops cannot be defaulted due to std::atomic<T> */
	TAnalysisProcess(TAnalysisProcess&& rhs) :
		q        {},
		_running ( false                  ), 
		_do_write( rhs._do_write          ),
		_thread  ( std::move(rhs._thread) ),
		_proc    ( std::move(rhs._proc)   ),
		info     ( std::move(rhs.info)    ),
		reader   ( std::move(rhs.reader)  ),
		writer   ( std::move(rhs.writer)  ) {}
	
	TAnalysisProcess& operator=(TAnalysisProcess&& rhs) noexcept {
		if(this == &rhs) return *this;
		_running  = false                 ;  
		_do_write = rhs._do_write         ; 
		_thread   = std::move(rhs._thread); 
		_proc     = std::move(rhs._proc)  ; 
		info      = std::move(rhs.info)   ; 
		reader    = std::move(rhs.reader) ; 
		writer    = std::move(rhs.writer) ; 
		return *this;
	}

	~TAnalysisProcess() = default;

	/* Adding a new proc mutates the type of the object. */

	/** 
	 * Reference-qualified: take an l/rvalue ref of an existing process object and
	 * eat up the object. Now we own it. Moves *this* object into a new one. 
	 */
	template <typename W,
		typename U = std::decay_t<W>
	> auto push_process(W&& w) && -> TAnalysisProcess<Ts..., U> {
		auto new_proc = std::tuple_cat(
			std::move(_proc),
			std::tuple<U>(std::move(w))
		);
		return TAnalysisProcess<Ts..., U> ( 
			std::move(new_proc), std::move(info), 
			std::move(reader), std::move(writer)
		);
	}
	
	/** 
	 * Emplace-style construction, again move *this* object into a new one. 
	 */
	template<typename U, typename... Args> 
	auto emplace_process(Args&&... args) && -> TAnalysisProcess<Ts..., U> {
		auto new_proc = std::tuple_cat(
			std::move(_proc),
			std::make_tuple( U(std::forward<Args>(args)... ) )
		);
		return TAnalysisProcess<Ts..., U> ( 
			std::move(new_proc), std::move(info), 
			std::move(reader), std::move(writer)
		);
	}

	/**
	 * A true independent copy of the object, meant to populate the pool singleton
	 * with clones of the initial processor object. 
	 */
	void Clone(TAnalysisProcess& dest) const { /* Only clone from the original object. */
		auto* p = std::get_if<std::unique_ptr<RExp::RNTupleParallelWriter>>(&writer.pwriter);
		/* ^^^ type: *std::unique_ptr<..> */
		if(!p) ERROR("Calling clone but original processor object is either unitialized or set to wrong state. "
				"State = %zu, 0 = Empty; 1 = Owning pointer; 2 = Raw pointer. Should be: " EMPH(1\n), writer.pwriter.index());

		dest._proc = this->_proc;
		dest.info  = this->info;
		dest.writer.pwriter = p->get();
		
		dest.Setup();
	}

	bool IsStopped() const noexcept { return ! _running && ! _thread.joinable(); }

	/**
	 * Returns the pointer 'T*' of the object wrapped up in the 'TAnalysisWorker<T>',
	 * at the tuple index 'I'
	 */
	template<std::size_t I>
	auto* GetProcess() noexcept {
		static_assert(I < Size(), "Request for process index outside of the tuple size.");
		return static_cast<base_t<I>*>(&std::get<I>(_proc));
	}
	template<std::size_t I>
	auto* GetProcess() const noexcept {
		static_assert(I < Size(), "Request for process index outside of the tuple size.");
		return static_cast<const base_t<I>*>(&std::get<I>(_proc));
	}
	// Runtime version, can throw.
	TProcessorBase* GetProcess(size_t i) {
		if(i >= Size()) ERROR("Request for worker index outside of the tuple size.");
		static constexpr auto table = make_getter_table(std::make_index_sequence<Size()>{});

		return table[i](this);
	}

	/**
	 * Return a table of processes, upcasted to the common parent.
	 */
	std::array<TProcessorBase*, Size()> GetProcesses() {
		std::array<TProcessorBase*, Size()> rv{};
		for(size_t i=0; i<Size(); ++i)
			rv[i] = this->GetProcess(i);
		return rv;
	}

	Int_t GetEntry(Long64_t entry) const { 
		if(util::IsEmpty(reader))
			ERROR("Empty input TTree/RNTuple. Invalid");
		else if(std::holds_alternative<util::RNPerThreadReader>(reader)) {
			std::get<util::RNPerThreadReader>(reader)
				._reader->LoadEntry(entry);
			return 0;
		}
		else {
			return std::get<util::TTreePerThreadReader>(reader)
			._tree->GetEntry(entry);
		}
	}

	u64 GetEntries() const { 
		if(util::IsEmpty(reader)) 
			ERROR("Empty input TTree/RNTuple in GetEntries call. Invalid");
		else if(std::holds_alternative<util::RNPerThreadReader>(reader)) {
			return static_cast<u64>( std::get <
				util::RNPerThreadReader
			> (reader)._reader->GetNEntries());
		}
		else {
			return static_cast<u64> (
				std::get<util::TTreePerThreadReader>(reader)._tree->GetEntries()
			);
		}
	}
	
	/* Run the setup. */
	void Setup() { SetupReader(); SetupWriter(); }

	void Start() {
		if(_running) 
			ERROR("Start called but worker thread is still marked as running? (%s)", _SELF_TYPE_CSTR);
		if(int v = writer.GetValidity(); v != 0) 
			ERROR("Calling start but writer handle isn't valid (v != 0). v = 0x%02x. Check API for GetValidity().", v); 
		if(int r = util::GetValidity(reader); r != 0)
			ERROR("Calling start but reader handle isn't valid (r != 0). r = 0x%02x. Check API for util::GetValidity().", r); 
		
		_running = true;
		
		WARN("Thread worker: TTree/RNTuple input: file %s, TFile* : 0x%016lx "
			"TTree* : 0x%016lx\n",
			info.in.fname.c_str(), (uintptr_t)std::get<util::TTreePerThreadReader>(reader)._file.get(),
			(uintptr_t)std::get<util::TTreePerThreadReader>(reader)._tree);
		std::apply([](auto&... ps) {
					(..., printf( "-- 0x%016lx [OUT:%s]\n-- 0x%016lx [IN:%s]\n", 
						(uintptr_t)ps.out.raw(), ps.out.GetName(), 
						(uintptr_t)std::get<0>(ps.in).raw(), std::get<0>(ps.in).GetName() ) );
				}, this->_proc);

		_thread = std::thread ( 
			[this] {
				util::Job j;
				while(_running.load(std::memory_order_relaxed)) {
					if( q.pop(j) ) {
						for(u64 evId = j.first; evId < j.last; ++evId) {
							this->GetEntry( static_cast<Long64_t>(evId) );
							
							std::apply([](auto&... ps) {
								(..., ps.ProcessEntry()); 
							}, this->_proc);
					
							if(this->_do_write)
								this->writer.ctx->Fill(*this->writer.entry);
						}
					} else {
#ifdef __HAS_SMALL_INTEL_SPIN
						_mm_pause(); /* Short pause, 100-150 clock cycles. */
#else
						{}           /* Do nothing; don't yield or reschedule - this is ~100 us latency. */
#endif
					}
				}

				/* On shutdown, drain the remainder of the queue. */
				while( q.pop(j) ) {
					for(u64 evId = j.first; evId < j.last; ++evId) {
						this->GetEntry( static_cast<Long64_t>(evId) );
						
						std::apply([](auto&... ps) {
								(..., ps.ProcessEntry()); 
						}, this->_proc);

						if(this->_do_write)
							this->writer.ctx->Fill(*this->writer.entry);	
					}
				}
			}
		);
	}

	void Stop() {
		_running.store(false, std::memory_order_relaxed);
		if(_thread.joinable()) _thread.join();
	}

	template<u32 N>
	auto MakePool(u32 NSlice) && -> TAnalysisPool<N, Ts...> {
		return TAnalysisPool<N, Ts...>( *this, NSlice );	
	}

	void Collect(const TAnalysisProcess& rhs) {
		std::tuple <
			std::pair<Ts&, const Ts&>...
		> pairs = util::zip_refs(this->_proc, rhs._proc);
		
		std::apply([](auto&... pr) {
				(..., pr.first.Collect( pr.second ));
			}, pairs
		);
	}
	
private:
	
	/* Each subthread writer is a slave to the initial one, who
	 * holds the true unique pointer handle. */
	void SetupWriter() {
		if(info.out.fname.empty() || info.out.out_rnname.empty())
			ERROR("Cannot proceed with setting up the writer if output (file,rnname) string information isn't given. (%s)", _SELF_TYPE_CSTR);
		if(writer.ctx || writer.entry) 
			ERROR("RNTupleWriter, context or entry already given (non-null)? (%s)", _SELF_TYPE_CSTR);

		RExp::RNTupleParallelWriter* pwriter_raw;
		
		/* Writer can be set up either from the original, which means we set up the bare model. */
		if(util::IsEmpty(writer.pwriter)) {
			auto model = RExp::RNTupleModel::CreateBare();
			
			util::for_each_in_tuple(this->_proc, [this, &model](auto& p /* TProcessor<Out(In...)>) */ ) 
				{
					const char* name = p.out.GetName();
					if(strlen(name) == 0) ERROR("Setting up the writer but an output container is unnamed. (%s)", _SELF_TYPE_CSTR);

					model->MakeField<typename decltype(p.out)::inner_type>( name );
				}
			);
			writer.pwriter = RExp::RNTupleParallelWriter::Recreate(std::move(model), info.out.out_rnname, info.out.fname);
			pwriter_raw = std::get<std::unique_ptr<RExp::RNTupleParallelWriter>> (writer.pwriter).get();
			if(!pwriter_raw)
				ERROR("RNTupleParallelWriter switched to original state, correct. But pointer is null after creation?");
		}
		else { /* ... Or it can be called from the clone. Initial `Clone()` call will set the raw pointer upfront. */
			if( ! std::holds_alternative<RExp::RNTupleParallelWriter*>(writer.pwriter) )
				ERROR("RNTupleParallelWriter isn't the original, but clone not switched to raw pointer handle. (%s)", _SELF_TYPE_CSTR);
			
			pwriter_raw = std::get<RExp::RNTupleParallelWriter*>(writer.pwriter);
			if(!pwriter_raw)
				ERROR("RNTupleParallelWriter switched to clone state, correct. But pointer is null?");
		}

		writer.ctx = pwriter_raw->CreateFillContext();
		writer.entry = writer.ctx->CreateEntry();

		util::for_each_in_tuple(this->_proc, [this](auto& p /* TProcessor<Out(In...)>) */ ) 
			{
				p.out._inner = this->writer.entry
					->GetPtr< typename decltype(p.out)::inner_type >( p.out.GetName() );
			}
		);
	}

	void SetupReader() {
		if(info.in.fname.empty())
			ERROR("Info given, but input file name empty. (%s)", _SELF_TYPE_CSTR);
		/* Based on the boolean inside, switch to either RN or TTree reader. */

		/* Reader object must be in empty state, */
		if(! util::IsEmpty(reader)) 
			ERROR("Trying to setup a fresh reader, but the object is already created with index: %zu. "
				"1 => RNReader; 2 => TTreeReader. (%s)", reader.index(), _SELF_TYPE_CSTR);
		
		std::string obj_name; bool is_ttree{};
		auto _file = std::make_unique<TFile>( info.in.fname.c_str() , "READ");
		if(!_file || _file->IsZombie() || !_file->IsOpen())
			ERROR("Setting up the reader but unable to make a TFile hook for \'%s\'. (%s)", info.in.fname.c_str(), _SELF_TYPE_CSTR);

		/* Try finding the read container. Can be either TTree or RNTuple. */
		for(TObject* _k : *_file->GetListOfKeys()) {
			TKey* k = dynamic_cast<TKey*>(_k);
			if(!k) continue;
			TClass* cl = TClass::GetClass(k->GetClassName());
			if(cl && cl->InheritsFrom(ROOT::RNTuple::Class()) ) {
				obj_name = std::string( k->GetName() ); is_ttree = false;
				break;
			} else if(cl &&  cl->InheritsFrom(TTree::Class()) ) {
				obj_name = std::string( k->GetName() ); is_ttree = true;
				break;
			}
		}

		if(obj_name.empty())
			ERROR("File %s opened fine, but cannot find a 'ROOT::RNTuple' or 'TTree' object inside? (%s)", 
				info.in.fname.c_str(), _SELF_TYPE_CSTR);
		
		_file.reset(nullptr);

		if(! is_ttree) {
			reader = util::RNPerThreadReader();
			auto& r = std::get<util::RNPerThreadReader>(reader);
			r._model = RExp::RNTupleModel::Create();
			
			/* For RNTupleModel' API, we use version 6.24.
			 * Newer ROOT 6.26 API could be different? */
			util::for_each_in_tuple(this->_proc, [this, &r](auto& p /* TProcessor<Out(In...)>) */ )
				{
					util::for_each_in_tuple(p.in, [this, &r](auto& cont /* In : TRawContainer */ )
						{
							/* Compile out the block if it's not a TContainer. */
							using ContType = typename std::remove_reference_t<decltype(cont)>;
							
							if constexpr( util::is_base_of_template<TContainer, ContType>::value ) {
								if(! strlen(cont.GetName()) )
									ERROR("Input container (RN-meant) unnamed? (%s)", _SELF_TYPE_CSTR);
								/* If the column has already been mapped, MakeField throws a 'RExp::RException'.
								 * In this case, just retrieve it and map it to the inner. */
								try {
									cont._inner = r._model->MakeField <
										typename ContType::inner_type
									> ( cont.GetName() );
								} catch(std::exception const& e) {
									cont._inner = r._model->GetDefaultEntry().GetPtr <
										typename ContType::inner_type
									> ( cont.GetName() );	
								} catch(...) {
									ERROR("Unknown exception caught? When trying to assing RNTuple column. (%s)",
										_SELF_TYPE_CSTR);
								}
							}
						}
					);
				}
			);

			r._reader = RExp::RNTupleReader::Open(std::move(r._model), info.in.fname, obj_name);
		}

		else { /* TTree version. */
			(void)TClass::GetClass("TTree");
			(void)TClass::GetClass("TChain");
			(void)TClass::GetClass("TBranch");
			(void)TClass::GetClass("TBasket");
			
			this->reader = util::TTreePerThreadReader();
			auto& r = std::get<util::TTreePerThreadReader>(reader);

			r._file = std::make_unique<TFile>( info.in.fname.c_str() , "READ");
			if(!r._file || r._file->IsZombie() || !r._file->IsOpen())
				ERROR("Unable to make a TFile hook for \'%s\'. (%s)", info.in.fname.c_str(), _SELF_TYPE_CSTR);

			r._tree = dynamic_cast<TTree*>(r._file->Get(obj_name.c_str()));
			if(!r._tree || r._tree->IsZombie())
				ERROR("Unable to make a TTree hook for \'%s\'. "
					"Even though TTree verified with name \'%s\'. (%s)", info.in.fname.c_str(), obj_name.c_str(), _SELF_TYPE_CSTR);
			
			r._tree->SetCacheSize(64*1024*1024);

			/* Ok, TTree API *sigh*.
			 * Namely, once we map a branch via `tree->SetBranchAddress(name, &ptr)`, then we can retrieve the pointers'
			 * address via: tree->GetBranch(name)->GetAddress . The return type is `char**` (???). 
			 * The only thing is that type safety isn't checked at runtime. Basically the other argument is `void**`
			 * somewhere down the line. This is really ugly, don't try this at home. 
			 * Checking the type does indeed work, if underlying type isn't templated. 
			 * It *should* demangle correctly. */

			WARN("~~~ New setup reader called, TFile* = 0x%016lx, TTree* = 0x%016lx <<<\n",
				(uintptr_t)r._file.get(), (uintptr_t)r._tree);
			util::for_each_in_tuple(this->_proc, [this, &r](auto& p /* TProcessor<Out(In...)>) */ )
				{
					util::for_each_in_tuple(p.in, [this, &r](auto& cont /* In : TRawContainer */ )
						{
							/* Compile out the block if it's not a TRawContainer. */
							using ContType = typename std::remove_reference_t<decltype(cont)>;

							if constexpr( util::is_base_of_template<TRawContainer, ContType>::value ) {
								if(! strlen(cont.GetName()) )
									ERROR("Input container (TTree-meant) unnamed? (%s)", _SELF_TYPE_CSTR);
								TBranch* b = r._tree->GetBranch( cont.GetName() );
								if(!b) 
									ERROR("File: \'%s\', TTree: \'%s\', branch \'%s\' not found. (%s)",
										this->info.in.fname.c_str(), r._tree->GetName(), cont.GetName(), _SELF_TYPE_CSTR);
								
								using T = typename ContType::inner_type;
								std::string c_type = util::type_name<T>();

								if(b->GetAddress()) { /* Means the branch is mapped already. */
									/* Try to check if types match. */
									const char* b_type = b->GetClassName();

									if(strcmp(b_type, c_type.c_str()) != 0)
										WARN("Setting up the container \'%s\'. "
											"File: \'%s\', TTree: \'%s\', branch \'%s\'. Types mismatch: "
											"TTree inspection gives branch type \'%s\', but we're requesting type \'%s\'. "
											"Is this correct? (%s)\n",
											cont.GetName(), this->info.in.fname.c_str(), r._tree->GetName(),
											b->GetName(), b_type, c_type.c_str(), _SELF_TYPE_CSTR);
									
									cont._inner = reinterpret_cast<T*>( *(void**)b->GetAddress() );
									WARN("Branch '%s', addr: 0x%016lx (TTree* addr: 0x%016lx), mapped to: 0x%016lx\n", 
										b->GetName(), (uintptr_t)b, (uintptr_t)r._tree, (uintptr_t)cont._inner); 
								} 
								else { // Branch isn't mapped yet. Map it now.
									/* Important to set to null. Else ROOT will not do anything, as it already
									 * dereferences to a valid T. 🤷 */
									cont._inner = nullptr; 
									
									(void)TClass::GetClass(c_type.c_str());
									
									r._tree->SetBranchAddress( cont.GetName(), &cont._inner);
									b->SetAutoDelete(kFALSE);
									WARN(">>>[N] '%s', addr: 0x%016lx (TTree* addr: 0x%016lx), mapped to: 0x%016lx\n", 
										b->GetName(), (uintptr_t)b, (uintptr_t)r._tree, (uintptr_t)cont._inner); 
								}
							} // if constexpr
						} 
					); // loop over input containers, per subprocess
				} 
			); // loop over subprocesses
		
			/* Warm-up. */
			Long64_t _n_entries = std::min(10LL, r._tree->GetEntries());
			for(Long64_t i=0; i<_n_entries; ++i) r._tree->GetEntry(i);

		} // TTree version end
	} // void SetupReader()
	
}; // TAnalysisProcess
