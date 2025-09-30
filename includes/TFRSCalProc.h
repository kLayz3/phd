#pragma once

#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TProcessor.h"

class TFRSCalProc : public TProcessor {
public:
	TFRSMapCont* input;
	TFRSCalCont* output;
	
	TFRSCalProc(TFRSMapCont& in, TFRSCalCont& out) : input(&in), output(&out) {}
	void ProcessEntry() noexcept;
	inline Int_t Write() override { return output->Write(); }

private:
	void ProcessTPC(int ) noexcept;
};
