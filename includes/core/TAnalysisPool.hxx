#pragma once

#include "AuxFunctions.hh"
#include "TROOT.h"
#include "libs.hh"
#include "TAnalysisProcess.hxx"
#include <thread>
#include <tuple>

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
#	warning "Running single-threaded.."
#	define POOL_MAX_THREADS 0
#endif

template <
	u32 N,
	typename... Processors
> struct TAnalysisPool final {
	static_assert(N <= POOL_MAX_THREADS_, 
		"TAnalysisProcess template instantiated with over-the-top capacity: " _TO_STRING(POOL_MAX_THREADS_) );
	static_assert(N >= 1, "TAnalysisProcess template parameter [1] size < 1? To run singlethreaded, "
		"put first integer template parameter to 1.");
	static_assert(!(N & (N-1)), "TAnalysisProcess template parameter [1] (n-processes) must be a power of two.");

	/* Only allowed ctor. */ 
	explicit TAnalysisPool(TAnalysisProcess<Processors...>& base, u32 _NSlice) : 
		pool{}, NSlice(_NSlice) 
	{
		assert((NSlice > 63) &&  "TAnalysisProcess constructor parameter [2] (slice size) must be bigger than 63 to be efficient.");
		pool[0] = std::move(base);
		pool[0].Setup();
		for(u32 i=1; i<N; ++i)
			pool[0].Clone(pool[i]);

		/* Once created, force reading of first few events to poke cling. */
		auto& w = this->Ref();
		w._do_write = false;
		
		u64 nLast = std::min (
			w.GetEntries(), (u64)10
		);
		
		volatile int sink = 0;
		for(u64 evId = 0; evId < nLast; ++evId) {
			w.GetEntry( static_cast<Long64_t>(evId) );
			sink += 1;
		}

		w._do_write = true;
	}

	/* On destructor sweep, write the single objects directly in the file. */
	~TAnalysisPool() { Collect(); Write(); g_loaded_containers.clear(); }

	/**
	 * Perform the dyadic fold of the Output containers, of each full process. 
	 * Process indexed with [0] contains the complete fold. Others are half-folded,
	 * and shouldn't be used any more. Idempotent function. 
	 */
	void Collect() {
		if(_is_collected) return;

		/* Check that all threads are safely merged. */
		for(const auto& p : pool)
			if(! p.IsStopped()) ERROR("Subthread isn't stopped while collector of the pool runs. Not allowed.");

		std::vector<
			TAnalysisProcess<Processors...> *
		> refs {};

		for(u32 i=0; i<N; ++i) refs.push_back( &pool[i] );
		dyadic_fold( refs );

		_is_collected = true;
	}
	
	/**
	 * Write the collected single- wise objects into the output TFile. Idempotent function. 
	 */
	void Write() {
		if( _is_written) return;

		/* Writing out the types for clarity. */
		TAnalysisProcess<Processors...>& process = Ref();
		util::IOInfo& info = process.info;
		
		/* First write the RNTuple. Gotta loop over all workers to close their
		 * RNTupleFillContext handles, only then destruct the RNTupleParallelWriter. */
		for(auto it = std::rbegin(pool); it != std::rend(pool); ++it) 
			it->writer.Reset();

		std::unique_ptr<TFile> f = std::make_unique<TFile>( info.out.fname.c_str(), "UPDATE");
		if(!f) 
			ERROR("Cannot open output file at end to write the single- wise objects. %s", info.out.fname.c_str());
		if(f->IsZombie()) 
			ERROR("Opened output file at end to write the single- wise objects is zombied?. %s", info.out.fname.c_str());
		if(!f->IsOpen()) 
			ERROR("Opened output file at end to write the single- wise objects, but is not `IsOpen()` %s", info.out.fname.c_str());
		if(!f->IsWritable())
			ERROR("Opened output file at end to write the single- wise objects, but is not writable? %s", info.out.fname.c_str());
		
		/* Fold the Write call over all subprocesses in pool[0] */
		util::for_each_in_tuple(process._proc /* tuple<TProcessor<Out(Ins..)> */, [&f](auto& subprocess)
			{
				/* subprocess: TProcessor<Out(Ins...)>& */
				for(auto& o : subprocess.out.GetTOnceVec()) {
					if(!o) ERROR("Requesting to write an object, but is nullptr? Subprocess type: %s",
						util::type_name<decltype(subprocess)>().c_str());
					o->Write( f.get() );
				}
			}
		);
		_is_written = true;
	}

	/**
	 * Send a batch of identical `NBatch` number of first entries 
	 * to each of the threads, to set up some initial parameters.
	 */
	void SendOneBatch(u64 startingIndex = 0, u32 NBatch = 0) {
		if(NBatch == 0) NBatch = NSlice;
		
		u64 nLast = std::min (
			pool[0].GetEntries(),
			(u64)NBatch + startingIndex
		);

		for(auto& w : pool) {
			w._do_write = false;
			w.Start();
		}

		for(u32 i=0; i<N; ++i) {
			int n_tries = 0;
			auto& w = pool[i];
			util::Job job { startingIndex, nLast };
			while(! w.q.push(job) ) {
				std::this_thread::yield();
				++n_tries;
			}
			/* Inspect and warn if somehow pushing a single job needed retries.. */
			if(n_tries != 0) WARN("SendOneBatch: needed %d retries for worker #%u\n", n_tries, i);
		}

		Stop(); /* This blocks until all the threads are joined. */
		for(auto& w : pool) 
			w._do_write = true;
	}

	void Start (
#ifdef __HAS_INDICATORS 
		indicators::ProgressBar& bar,
#endif
		const u64 max_entries = static_cast<u64>(-1)
	) {
		auto& ref_process = Ref();

		u64 nentries = std::min (
			ref_process.GetEntries(),
			max_entries
		);

		ROOT::EnableThreadSafety();
		for(auto& w : pool)
			w.Start();

#ifdef __HAS_INDICATORS
		indicators::show_console_cursor(false);
#endif
		u32 next = 0;
		for(u64 i = 0; i < nentries; i+=NSlice) {
			util::Job j {
				i,                             // Starting index, included.
				std::min(i + NSlice, nentries) // Stopping index, excluded.
			};

			/* Choose a process thread; round-robin. 
			 * If currently selected thread has full queue, do a short spin and try next one. */
			for(auto& w = pool[next]; ! w.q.push(j); ++next, next %= N) {
#ifdef __HAS_SMALL_INTEL_SPIN
				_mm_pause(); /* Short pause, 100-150 clock cycles. */
#else
				{}           /* Do nothing; don't yield or reschedule - this is ~100 us latency. */
#endif
			}
#ifdef __HAS_INDICATORS
			util::PrintProgress(bar, j.last, nentries, NSlice-1 );
#endif
		}

		Stop();

#ifdef __HAS_INDICATORS
		bar.mark_as_completed();
		indicators::show_console_cursor(true);
#endif

	} // void Start(...)

	void Stop() {
		for(auto& w : pool) w.Stop(); 
	}
	
	decltype(auto) GetPool()       noexcept { return ( this->pool ); }
	decltype(auto) GetPool() const noexcept { return ( this->pool ); }
	
	TAnalysisProcess<Processors...>&       Ref()       noexcept { return pool[0]; }
	TAnalysisProcess<Processors...> const& Ref() const noexcept { return pool[0]; }

private:
	std::array <
		TAnalysisProcess<Processors...>, N
	> pool;
	u32 NSlice;

	bool _is_collected { false };
	bool _is_written   { false };

	template<typename T>
	void dyadic_fold(std::vector<T*> v) {
		size_t Nv = v.size();
		if(Nv & (Nv-1)) ERROR("Dyadic fold container size ill-formed, is %zu, but should be power of 2.", Nv);

		if(Nv == 1) return;
		
		const size_t half = Nv / 2;
		std::vector<T*> next(half);

		for(size_t i=0; i<half; ++i) {
			v[ 2*i ] -> Collect( (const T&)(*v[2*i + 1]) );
			next[i] = v[2*i];
		}

		dyadic_fold(next);
	}
};

/* ===============================================================
 * =============================================================== */
/* Invariant API for both multithreaded and singlethreaded modes.  */
/* NSlice field is not really used, but still keep it to keep identical calls across the
 * specialization. */

template<typename... Processors>
struct TAnalysisPool<1, Processors...> final {
	explicit TAnalysisPool(TAnalysisProcess<Processors...>& base, u32 _NSlice) : 
		pool{}, NSlice(_NSlice) 
	{
		pool[0] = std::move(base);
		pool[0].Setup();
	}

	/* On destructor sweep, write the single objects directly in the file. */
	~TAnalysisPool() { Collect(); Write(); }
	
	void Collect() {}

	/**
	 * Write the collected single- wise objects into the output TFile. Idempotent function. 
	 */
	void Write() {
		if( _is_written) return;

		/* Writing out the types for clarity. */
		TAnalysisProcess<Processors...>& process = Ref();
		util::IOInfo& info = process.info;
		
		/* First write the RNTuple. */
		process.writer.Reset();	

		std::unique_ptr<TFile> f = std::make_unique<TFile>( info.out.fname.c_str(), "UPDATE");
		if(!f || f->IsZombie()) 
			ERROR("Cannot open output file at end to write the single- wise objects. %s", info.out.fname.c_str());
		if(!f->IsOpen()) 
			ERROR("Opened output file at end to write the single- wise objects, but is not `IsOpen()` %s", info.out.fname.c_str());
		if(!f->IsWritable())
			ERROR("Opened output file at end to write the single- wise objects, but is not writable %s", info.out.fname.c_str());
		
		/* Fold the Write call over all subprocesses in pool[0] */
		util::for_each_in_tuple(process._proc /* tuple<TProcessor<Out(Ins..)> */, [&f](auto& subprocess)
			{
				/* subprocess: TProcessor<Out(Ins...)>& */
				for(auto& o : subprocess.out.GetTOnceVec()) {
					if(!o) ERROR("Requesting to write an object, but is nullptr? Subprocess type: %s",
						util::type_name<decltype(subprocess)>().c_str());
					o->Write( f.get() );
				}
			}
		);
		_is_written = true;
	}
	
	/**
	 * Send a batch of identical `NBatch` number of first entries 
	 * to the underlying process, to set up some initial parameters.
	 */
	void SendOneBatch(u64 startingIndex = 0, u32 NBatch = 0) {
		if(NBatch == 0) NBatch = (NSlice > 0) ? NSlice : 2048;
		auto& process = Ref();
		u64 nLast = std::min (
			process.GetEntries(),
			(u64)NBatch + startingIndex
		);
		WARN("A quick singlethreaded, sending a batch and processing! NSubProc = %zu\n", process.Size());	
		for(u64 evId = startingIndex; evId < nLast; ++evId) {
			process.GetEntry( static_cast<Long64_t>(evId) );

			std::apply([](auto&... ps) {
					(..., ps.ProcessEntry()); 
				}, process._proc
			);
			/* Don't fill. */
		}
	}

	/* In singlethreaded case, don't poke the underlying std::thread of the process,
	 * just call ProcessEntry() directly. */
	void Start (
#ifdef __HAS_INDICATORS 
		indicators::ProgressBar& bar,
#endif
		const u64 max_entries = static_cast<u64>(-1)
	) {
		auto& process = Ref();

		u64 nentries = std::min (
			process.GetEntries(),
			max_entries
		);
		u64 n_print_every = ((NSlice > 0) ? NSlice : 512); 
		for(u64 evId = 0; evId < nentries; ++evId) {
			process.GetEntry( static_cast<Long64_t>(evId) );

#ifdef __HAS_INDICATORS
			util::PrintProgress(bar, evId, nentries, n_print_every);
#endif
			std::apply([](auto&... ps) {
					(..., ps.ProcessEntry()); 
				}, process._proc
			);
			process.writer.ctx->Fill( *process.writer.entry );
		}
	}
	
	void Stop() { Ref().Stop(); } /* No-op; define it anyway to keep API invariant. */

	decltype(auto) GetPool()       { return ( this->pool ); }
	decltype(auto) GetPool() const { return ( this->pool ); }

	TAnalysisProcess<Processors...>&       Ref()       noexcept { return pool[0]; }
	TAnalysisProcess<Processors...> const& Ref() const noexcept { return pool[0]; }
	
private:
	std::array <
		TAnalysisProcess<Processors...>, 1
	> pool;
	u32 NSlice;

	/* bool _is_collected { false }; */
	bool _is_written   { false };
};
