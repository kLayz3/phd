#include "TFRSMapCont.h"
#include "Rtypes.h"

TFRSMapCont::TFRSMapCont() : TContainer("FRS") {}


void TFRSMapCont::Setup() {
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
			Form("TPC%d_MR", i), Form("Multipl. TPC%d right dl(0)", i), 10, 0, 10
		);
		h1_tpc_ma1[i] = RegisterObject<TH1I>(
			Form("TPC%d_MA1", i), Form("Multipl. TPC%d anode(1)", i), 10, 0, 10
		);
		h1_tpc_ma2[i] = RegisterObject<TH1I>(
			Form("TPC%d_MA2", i), Form("Multipl. TPC%d anode(2)", i), 10, 0, 10
		);
		for(int a=0; a<4; ++a) {
			h1_tpc_csum[i][a] = RegisterObject<TH1I>(
				Form("TPC%d_CSum%d", i,a), Form("CSum TPC%d anode(%d) mult == 1", i, a), 20000, 10000,30000
			);
			h1_tpc_ydiff[i][a] = RegisterObject<TH1I>(
				Form("TPC%d_YDiff%d", i,a), Form("TPC%d anode(%d) - ref (mult == 1; in both)", i, a), 20000, 0, 100000
			); 
		}
	}
};

ClassImp(RNSciMap);
ClassImp(RNTPCMap);
ClassImp(RNTPCMap::Measurement);
ClassImp(RNMUSICMap<8>);
ClassImp(RNFRSMap);
