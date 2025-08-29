#pragma once

/* Main class used for FOOT analysis and removing the pedestals. NOTE: this class is ONLY used for offline,
 * as it makes huge performance overheads to make a separate entry loop per each instance of the class.
 * This was done to help the code readability when exported to different other projects. */

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TTree.h"
#include "libs.hh"
#include "TProcessor.h"
#include "TFOOTPedestalCont.h"
#include "TOnce.hxx"
#include <utility>

class TFOOTPedestalProc : public TProcessor {
public:
	static constexpr int N_STRIPS          = TFOOTPedestalCont::N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = TFOOTPedestalCont::N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = TFOOTPedestalCont::N_STRIPS_PER_ASIC; /* 64  */
	
	static constexpr double BAD_STRIP_CUTOFF = 3.0;	
public:
	enum ProcessType {
		kGPED,
		kEPED
	} process_type = kGPED;
		
	TFOOTPedestalCont* data;

	TFOOTPedestalProc(TFOOTPedestalCont& data) :  
		data(&data) {}

	/* Rule-of-five: if a destructor or a custom move ctor/assignment
	 * is declared, then all three (dtor/move ctor/move assignment) must be declared,
	 * even if `= default;`. The copies are implicitly deleted, due to the base class. */
	 /* 
	TFOOTPedestalProc(TFOOTPedestalProc&& ) noexcept = default;
	TFOOTPedestalProc& operator=(TFOOTPedestalProc&& ) noexcept = default;
	~TFOOTPedestalProc() = default;
	*/

	void ProcessGlobalPedestal() noexcept;
	void CalcGlobalPedestal();
	void ProcessEventPedestal() noexcept;
	void CalcFinalPedestal();

	void ProcessEntry() noexcept ;

	Int_t Write() override { return data->Write(); }

	int GetRawADC(int chn);
	static constexpr inline int GetChannel(int asic, int offset) { return asic * N_STRIPS_PER_ASIC + offset; }

private:
	std::array<double, N_STRIPS_PER_ASIC> ped_offset = {0.0};
};
