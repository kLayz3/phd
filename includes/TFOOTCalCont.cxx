#include "TFOOTCalCont.h"
#include "Rtypes.h"
#include "RtypesCore.h"
#include "TH1D.h"

TFOOTCalCont::TFOOTCalCont() : FOOT_N(-1), POS(-1), N(0) {
	fCX.reserve(INIT_CAPACITY);
	fCE.reserve(INIT_CAPACITY);
	fCM.reserve(INIT_CAPACITY);
	fCT.reserve(INIT_CAPACITY);
}

TFOOTCalCont::~TFOOTCalCont() {}

void TFOOTCalCont::Init(TDictInfo info) {
	std::unordered_map<const char*, int> int_mappings;
	for(auto key : {"FOOT_ID", "FOOT_POS"}) {
		auto n_it = info.find(key);
		if(n_it == info.end())
			ERROR("\'%s\' key not found in the info hashmap.", key);
		int n;
		try {
			n = std::stoi(n_it->second);
		} catch(const std::exception& e) {
			ERROR("For key \'%s\', value found: " EMPH(%s) " but unparsable to integer. Err: %s", key, n_it->second.c_str(), e.what());
		}
		int_mappings.emplace(key, n);
	}

	FOOT_N = int_mappings.at("FOOT_ID");
	POS = int_mappings.at("FOOT_POS");

	this->SetName(Form("CFOOT%d", POS));

	assert(GetOwnedTOnceObjects().size() == 0 && "Don't call `SetId` twice. Clear the owned objects first.");
	
	h1_raw_mult = RegisterObject<TH1I>("h1_raw_mult", Form("FOOT(%2d: %d) raw multiplicity", FOOT_N, POS), 50,0,50);
	h1_mult     = RegisterObject<TH1I>("h1_mult", Form("FOOT (%2d: %d) multiplicity", FOOT_N, POS), 200,0,50);
	h1_dE       = RegisterObject<TH1I>("h1_dE", Form("FOOT(%2d: %d) energy", FOOT_N, POS), 5000,0,1000);
	h1_X        = RegisterObject<TH1I>("h1_x", Form("FOOT(%2d: %d) clust position", FOOT_N, POS), 5*N_STRIPS, 0, N_STRIPS);
	h1_cl_type  = RegisterObject<TH1I>("h1_cl_type", Form("(%2d:%d) clust type", FOOT_N, POS), 5, 0, 5);
	h1_dE_m1    = RegisterObject<TH1I>("h1_dE_m1", Form("(%2d:%d) energy with multiplicity 1", FOOT_N, POS), 2500, 0, 500);
	h1_dE_m2    = RegisterObject<TH1I>("h1_dE_m2", Form("(%2d:%d) energy with multiplicity 2", FOOT_N, POS), 2500, 0, 500);
	h1_dE_m3    = RegisterObject<TH1I>("h1_dE_m3", Form("(%2d:%d) energy with multiplicity 3", FOOT_N, POS), 2500, 0, 500);

	_fBadE.reserve(N_STRIPS);
}

void TFOOTCalCont::Clean(Option_t* option) noexcept {
	(void)option;
	N = 0;
	fCX.clear();
	fCE.clear();
	fCM.clear();
	fCT.clear();
	_fBadE.clear();
}

void TFOOTCalCont::AddCluster(double mean_x, double e, double mult, ClusterType ct) {
	fCX.push_back(mean_x);
	fCE.push_back(e);
	fCM.push_back(mult);
	fCT.push_back(ct);
	++N;
}

IMPL_CONTAINER_METHODS(TFOOTCalCont)

ClassImp(TFOOTCalCont);
