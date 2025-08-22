#pragma once

#include <mutex>
#include <atomic>
#include <thread>
#include <condition_variable>
#include "TNamed.h"
#include "TContainer.h"
#include "libs.hh"
#include "TOnceBase.h"

/* This is a base class of analysis workers.
 * Single TAnalysisPool has a handle for multiple of them, and will call their virtual DTOR's
 * properly at the end. */
/* Users must extend this class, else it's a no-op. Own classes will either hold data directly or via refs' & ptrs'.
 * Either way, the underlying worker thread _thread will have all the available handles via `this` ptr. */

class alignas(64) TAnalysisWorker {
	template<size_t N> friend class TAnalysisPool;
	std::thread _thread; //!
private:
	std::atomic<bool> _stop;     //!
	std::atomic<bool> _has_work; //!
public:
	TContainer& input;
	TContainer& output;
	
	TAnalysisWorker(TContainer& in, TContainer& out) : input(in), output(out) {}
	TAnalysisWorker() = delete;
	TAnalysisWorker(const TAnalysisWorker& ) = delete;
	TAnalysisWorker(TAnalysisWorker&& ) = delete;
	TAnalysisWorker& operator=(const TAnalysisWorker& ) = delete;
	TAnalysisWorker& operator=(TAnalysisWorker&& ) = delete;

	virtual ~TAnalysisWorker(); 

	inline bool IsDone() const noexcept { return ! _has_work.load(std::memory_order_acquire); }
	
	/**
	 * Users must override these method, else it's a no-op. 
	 */
	virtual void Init(const TDictInfo& , const TDictInfo& );
	virtual void ProcessEntry() noexcept;
	virtual void Clear(Option_t* option = "") noexcept;
	
	void Start();

public: // related to parts to be serialized only once.
	
	/**
	 * Serialize single objects (non-TTree) parts of the worker to the currently open directory.
	 */
	virtual Int_t Write();
};
