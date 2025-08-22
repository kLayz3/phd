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
#include "TAnalysisWorker.h"
#include "TFOOTPedestalCont.h"
#include "TOnce.hxx"
#include <utility>

class TFOOTPedestalProc : public TAnalysisWorker {
public:
	static constexpr int N_STRIPS          = TFOOTPedestalCont::N_STRIPS; //!
	static constexpr int N_ASIC            = TFOOTPedestalCont::N_ASIC; //!
	static constexpr int N_STRIPS_PER_ASIC = TFOOTPedestalCont::N_STRIPS_PER_ASIC; //!

public:
	enum ProcessType {
		kGPED,
		kEPED
	} process_type = kGPED;
		
	TFOOTPedestalCont& data;
	TFOOTPedestalProc(TFOOTPedestalCont& data) : TAnalysisWorker(data, data), 
		data(data) {}
	~TFOOTPedestalProc() = default;

	void ProcessGlobalPedestal() noexcept;
	void CalcGlobalPedestal();
	void ProcessEventPedestal() noexcept;
	void CalcFinalPedestal();

	void ProcessEntry() noexcept override;

	int GetRawADC(int chn);
	static constexpr inline int GetChannel(int asic, int offset) { return asic * N_STRIPS_PER_ASIC + offset; }

	inline void Clear(Option_t* option = "") noexcept override { data.Clean(option); }

private:
	std::array<double, N_STRIPS_PER_ASIC> ped_offset = {0.0};
};
