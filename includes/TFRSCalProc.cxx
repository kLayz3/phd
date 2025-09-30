#include "TFRSCalProc.h"

using TPCParam = TFRSCalCont::TPCParam;

void TFRSCalProc::ProcessEntry() noexcept {
	for(int i=0; i<7; ++i)
		this->ProcessTPC(i);				
}

void TFRSCalProc::ProcessTPC(int i) noexcept {	
	TPCParam const& p = (*output->tpc_param)[i];
	RNTPCMap const& in = input->inner().tpc[i];
} 
