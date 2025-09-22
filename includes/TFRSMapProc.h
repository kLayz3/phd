#pragma once

#include "TProcessor.h"
#include "TFRSMapCont.h"

class TFRSMapProc : public TProcessor {
public:
	TFRSMapCont* data;
	TFRSMapProc(TFRSMapCont& data, int do_analysis = 1) : 
		data(&data),
		do_analysis(do_analysis) {}

	void  ProcessEntry() noexcept; 
	void _ProcessEntry() noexcept;
	
	int do_analysis = 1; /* 1 or 0 */
	inline Int_t Write() override { return data->Write(); }
};
