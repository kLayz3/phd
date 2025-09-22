#include "TFRSMapCont.h"

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

TFRSMapCont::~TFRSMapCont() {}

u32 TFRSMapCont::TPC::tpc_ref_i[]{
	0, // TPC21 => SCI21
	0, // TPC22 => SCI21
	0, // TPC23 => SCI21
	1, // TPC24 => SCI22
	2, // TPC31 => SCI31
	3, // TPC41 => SCI41
	3  // TPC42 => SCI41
};

void TFRSMapCont::Clear() { 
	for(int i=0; i<4; ++i)
		sci[i].__clear();

	for(int i=0; i<7; ++i)
		tpc[i].__clear();

	for(int i=0; i<2; ++i)
		music[i].__clear();

	tpat = 0;
}

IMPL_CONTAINER_METHODS(TFRSMapCont)

ClassImp(TFRSMapCont::Sci)
ClassImp(TFRSMapCont::TPC)
ClassImp(TFRSMapCont::MUSIC<8>)
ClassImp(TFRSMapCont)
