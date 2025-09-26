#include "TFRSMapCont.h"
#include "Rtypes.h"

TFRSMapCont::TFRSMapCont() : TContainer("FRS") {
	for(int i=0; i<4; ++i) {
		h1_sci_ml[i] = RegisterObject<TH1I>(
			Form("SCI%d_ML", i), Form("Multipl. SCI%d left", i), 10, 0, 10
		);
		h1_sci_mr[i] = RegisterObject<TH1I>(
			Form("SCI%d_MR", i), Form("Multipl. SCI%d right", i), 10, 0, 10
		);
	}

	for(int i=0; i<7; ++i) {
		h1_tpc_ml[i] = RegisterObject<TH1I>(
			Form("TPC%d_ML", i), Form("Multipl. TPC%d left dl(0)", i), 10, 0, 10
		);
		h1_tpc_mr[i] = RegisterObject<TH1I>(
			Form("TPC%d_ML", i), Form("Multipl. TPC%d right dl(0)", i), 10, 0, 10
		);
		h1_tpc_ma1[i] = RegisterObject<TH1I>(
			Form("TPC%d_MA1", i), Form("Multipl. TPC%d anode(1)", i), 10, 0, 10
		);
		h1_tpc_ma2[i] = RegisterObject<TH1I>(
			Form("TPC%d_MA2", i), Form("Multipl. TPC%d anode(2)", i), 10, 0, 10
		);
	}
};

ClassImp(RNSciMap);
ClassImp(RNTPCMap);
ClassImp(RNMUSICMap<8>);
ClassImp(RNFRSCont);
