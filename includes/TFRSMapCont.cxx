#include "TFRSMapCont.h"

TFRSMapCont::TFRSMapCont() : TContainer("FRS") {}

void TFRSMapCont::Setup() {
	for(int i=0; i<N_VALID_SCI; ++i) {
		h1_sci_ml[i] = RegisterObject<TH1I>(
			Form("SCI%d_ML", i), Form("Multipl. SCI%d left", i), 10, 0, 10
		);
		h1_sci_mr[i] = RegisterObject<TH1I>(
			Form("SCI%d_MR", i), Form("Multipl. SCI%d right", i), 10, 0, 10
		);
		h1_sci_diff_lr[i] = RegisterObject<TH1I>(
			Form("SCI%d_diff_lr", i), Form("TDC diff. SCI%d L-R (multp=1)", i), 2000, -1000, 1000
		);
	}

	for(int i=0; i<N_VALID_TPC; ++i) {
		for(int d=0; d<2; ++d) {
			h1_tpc_dl_lraw[i][d] = RegisterObject<TH1I>(
				Form("TPC%d_DL%d", i,d), Form("TPC%d delayline(%d) left raw", i, d), 20000, 0, 60000
			);
			h1_tpc_dl_rraw[i][d] = RegisterObject<TH1I>(
				Form("TPC%d_DR%d", i,d), Form("TPC%d delayline(%d) right raw", i, d), 20000, 0, 60000
			);
		}
		for(int a=0; a<4; ++a) {
			h1_tpc_araw[i][a] = RegisterObject<TH1I>(
				Form("TPC%d_A%d", i, a), Form("TPC%d anode(%d) raw", i, a), 20000, 0, 60000
			);
		}
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

		for(int a=0; a<2; ++a) {
			h1_tpc_adiff[i][a] = RegisterObject<TH1I>(
				Form("TPC%d_ADiff%d", i,a), Form("TPC%d anode(%d) - anode(%d) (mult == 1; in both)", i, 2*a, 2*a+1), 
					20000, -10000, 10000
			); 
		}
		h1_tpc_ldiff[i] = RegisterObject<TH1I>(
			Form("TPC%d_LDiff", i), Form("TPC%d DL0 left - DL1 left (mult == 1; in both)", i), 20000, -10000, 10000
		); 
		h1_tpc_rdiff[i] = RegisterObject<TH1I>(
			Form("TPC%d_RDiff", i), Form("TPC%d DL0 right - DL1 right (mult == 1; in both)", i), 20000, -10000, 10000
		); 
	}
};

TTrigMapCont::TTrigMapCont() : TContainer("Trig") {}

void TTrigMapCont::Setup() {
	h1_full_tpat = RegisterObject<TH1I>("full_tpat", "Complete Trigger Pattern", 65536, 0, 65536);
	h1_tpat = RegisterObject<TH1I>("tpat", "Trigger Pattern", 17, 0, 17);
	h1_wr_diff = RegisterObject<TH1I>("wr_diff", "FRS Whiterabbit Increment [us]", 200, -5, 2000);
}

ClassImp(RNSciMap);
ClassImp(RNSciMap::Measurement);
ClassImp(RNTPCMap);
ClassImp(RNTPCMap::Measurement);
ClassImp(RNMUSICMap<8>);
ClassImp(RNFRSMap);

ClassImp(RNTrigMap);
