#pragma once

#include "libs.hh"
#include "dbg.hh"
#include "TTree.h"
#include "TChain.h"
#include <atomic>
#include <variant>
#include <optional>
#include "TAnalysisWorker.h"
#include "TContainer.h"

/* Passed from build tool. */
#ifndef POOL_MAX_THREADS_
#	define POOL_MAX_THREADS_ 10
#endif

class TAnalysisWorker;

template<size_t MAX_WORKERS>
class TAnalysisPool final {
	static_assert(MAX_WORKERS <= POOL_MAX_THREADS_, 
		"TAnalysisPool template instantiated with over-the-top capacity. Consider lowering.");
public:
	constexpr size_t Size() { return MAX_WORKERS; }

private:
	std::array<TAnalysisWorker*, MAX_WORKERS> _pool = {nullptr}; //!
	int n_valid_workers = 0; //!

private:
	bool _stop = true; //! Flagging this will join the workers' threads back to the main
	std::array<bool, MAX_WORKERS> _owned = {0}; //! A label indicator which Worker object should be owned by the pool.	

public:
	// TTree objects are owned and handled by gROOT. Here we expose just the raw handles.
	std::variant<TTree*, TChain*> in; //!
	std::optional<TTree*> out; //!

	TAnalysisPool() = default;
	TAnalysisPool(const TAnalysisPool& ) = delete;
	TAnalysisPool(TAnalysisPool&& ) = delete;
	TAnalysisPool& operator=(const TAnalysisPool& ) = delete;
	TAnalysisPool& operator=(TAnalysisPool&& ) = delete;
	~TAnalysisPool();
	
	inline bool IsStopped() const noexcept { return _stop; }
	inline TAnalysisWorker* GetWorker(size_t i) const noexcept {
		if(i >= MAX_WORKERS) return nullptr;
		else return _pool[i];
	}
	inline Int_t GetEntry(Long64_t entry, Int_t getall = 0) const { 
		if(in.valueless_by_exception()) 
			ERROR("Valueless input TTree/TChain. Invalid");
		else return std::visit([=](const auto& obj) { return obj->GetEntry(entry,getall); }, this->in); 
	}
	inline Long64_t GetEntries() const noexcept { 
		if(std::holds_alternative<TTree*>(in)) 
			return std::get<TTree*>(in)->GetEntries();
		else if(std::holds_alternative<TChain*>(in)) 
			return std::get<TChain*>(in)->GetEntriesFast();
		else return 0;
	}
	inline Int_t FillOutput() const noexcept { 
		if(out and *out) return (*out)->Fill();
		else return 0;
	}
	
	/**
	 * Add a base class ptr of a either stack or heap-allocated worker, to be managed by the pool.
	 * Every action the worker needs to take is completely managed by the pool instance.
	 * The pool however will not free the object. 
	 */
	TAnalysisPool& AddWorker(TAnalysisWorker*);

	/**
	 * Add a heap-allocated worker to be owned by the pool. It will call destructors on release.
	 */
	TAnalysisPool& AddOwnedWorker(TAnalysisWorker*);

	/**
	 * Start the underlying thread of every worker.
	 */
	void Start();

	/**
	 * Mark the flag to wake-up for all workers simultaneously.
	 */
	void AssignWork() const noexcept;

	/**
	 * Fill the remainder of the pool with no-op workers. 
	 */
	void FinalizeInit();

	/**
	 * Block the master thread until the execution of all the workers is marked finished. 
	 */
	void Await() const noexcept;

	/**
	 * Stops execution of all the workers and joins their internal threads.
	 */
	void Stop();

	/**
	 * Stops execution and clears all the workers. Call the underlying destructor of owned workers.
	 */
	void ClearPool();
};

/* Explicitly check for nullptr here to allow this call tailing `ClearPool`. */
template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::Stop() {
	_stop = true;
	for(size_t i=0; i < MAX_WORKERS; ++i) {
		if(_pool[i] != nullptr and _pool[i]->_thread.joinable()) {
			_pool[i]->_stop.store(true, std::memory_order_release);
			_pool[i]->_thread.join();
		}
	};
}

template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::ClearPool() {
	this->Stop();

	for(size_t i=0; i < MAX_WORKERS; ++i) {
		if(_owned[i]) delete _pool[i];
		_pool[i] = nullptr;
	}
	n_valid_workers = 0;
}

template<size_t MAX_WORKERS>
TAnalysisPool<MAX_WORKERS>::~TAnalysisPool() {
	this->ClearPool();
}

template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::Await() const noexcept {
	while(true) {
		bool all_done = 1;

		for(size_t i=0; i < MAX_WORKERS; ++i)
			all_done &= _pool[i]->IsDone();
		
		if(all_done) break; 
	}
}

template<size_t MAX_WORKERS>
TAnalysisPool<MAX_WORKERS>& TAnalysisPool<MAX_WORKERS>::AddWorker(TAnalysisWorker* w) {
	if(n_valid_workers >= static_cast<int>(MAX_WORKERS))
		ERROR("Failed assert (n_valid_workers < MAX_WORKERS): (%d<%d)\n", n_valid_workers, static_cast<int>(MAX_WORKERS));
	
	if(!w) ERROR("Worker arg (0x%lx) must not be nullptr.\n", (uintptr_t)w);
	_pool[n_valid_workers++] = w;
	return *this;
}

template<size_t MAX_WORKERS>
TAnalysisPool<MAX_WORKERS>& TAnalysisPool<MAX_WORKERS>::AddOwnedWorker(TAnalysisWorker* w) {
	this->AddWorker(w);
	_owned[n_valid_workers - 1] = true;
	return *this;
}

template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::AssignWork() const noexcept {
	for(size_t i=0; i<MAX_WORKERS; ++i)
		_pool[i]->_has_work.store(true, std::memory_order_release);
}

template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::FinalizeInit() {
	while(n_valid_workers < static_cast<int>(MAX_WORKERS)) {
		this->AddOwnedWorker(new TAnalysisWorker(TContainer::dummy, TContainer::dummy));
		WARN("Added a no-op worker: %i; consider lowering the capacity of TAnalysisPool instance (%zu).\n", n_valid_workers, MAX_WORKERS);
	}
}


template<size_t MAX_WORKERS>
void TAnalysisPool<MAX_WORKERS>::Start() {
	/* State prior to this call must be stopped. */
	if(!_stop) ERROR("Must be in a stopped state before calling start.");
	_stop = false;

	for(size_t i=0; i<MAX_WORKERS; ++i) {
		if(!_pool[i]) ERROR("Worker[%zu] (0x%lx) must not be nullptr.\n", i, (uintptr_t)_pool[i]);
		_pool[i]->_stop.store(false, std::memory_order_release);
		_pool[i]->Start();
	}
}

