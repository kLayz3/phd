#include "TAnalysisWorker.h"
#include "AuxFunctions.hh"
#include <atomic>

TAnalysisWorker::~TAnalysisWorker() {
	if(_thread.joinable()) _thread.join();
}
// Users can optionally override these methods.

// Just call the init of underlying ref. Input should be just tied to inputting from previous step.
void TAnalysisWorker::Init(const TDictInfo& info_in, const TDictInfo& info_out) { 
	input.Init(info_in);
	output.Init(info_out);
}
void TAnalysisWorker::Clear(Option_t* option) noexcept { output.Clean(); }

// This method should be overwritten. 
void TAnalysisWorker::ProcessEntry() noexcept {}

void TAnalysisWorker::Start() {
	if(_thread.joinable())
		ERROR("(%s => %s) trying to (re)start the worker but thread is in a joinable state?", input.GetName(), output.GetName());
	if(_stop.load(std::memory_order_acquire)) 
		ERROR("(%s => %s) trying to start the worker but it's marked as stopped?", input.GetName(), output.GetName());

	_thread = std::thread{
		[this]() {
			while(! _stop.load(std::memory_order_acquire)) {
				if(_has_work.load(std::memory_order_acquire)) {
					ProcessEntry();
					_has_work.store(false, std::memory_order_release);
				}
			}
		}
	};
}

Int_t TAnalysisWorker::Write() {
	return output.Write();
}
