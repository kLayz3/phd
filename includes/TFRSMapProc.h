#pragma once

#include "TProcessorBase.h"
#include "TFRSMapCont.h"

struct TFRSMapProc : TProcessorBase {
public:
	TFRSMapCont* data;
	TFRSMapProc(TFRSMapCont& data, int do_analysis = 1) : 
		data(&data),
		do_analysis(do_analysis) {}

	void  ProcessEntry() noexcept; 
	void _ProcessEntry() noexcept;
	
	int do_analysis = 1; /* 1 or 0 */
	inline Int_t Write() override { return data->Write(); }
private:
	/* Some TPC vars' local temporary storage. */
	Int_t _nhits_l[2] {};
	Int_t _nhits_r[2] {};
	Int_t _nhits_a[4] {};
	Int_t _nhits_s {};
	std::array<Int_t, RNTPCMap::MAX_SIZE> _temp {};
};
