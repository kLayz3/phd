#include "TFOOTCalCont.h"
#include "Rtypes.h"
#include "RtypesCore.h"
#include "TH1D.h"

TFOOTCalCont::TFOOTCalCont() {
	fCX = &_x[0];
	fCE = &_e[0];
	fCM = &_m[0];
	fCT = &_t[0]; 
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

	this->SetName(Form("FOOT%d", POS));

	assert(GetOwnedTOnceObjects().size() == 0 && "Don't call `SetId` twice. Clear the owned objects first.");
	
	h1_mult = RegisterObject<TH1D>("h1_mult", Form("(%d -> %d) multiplicity", FOOT_N, POS), 100,0,20);
	h1_dE   = RegisterObject<TH1D>("h1_dE", Form("(%d -> %d) multiplicity", FOOT_N, POS), 100,0,20);
	h2_X    = RegisterObject<TH1D>("h2_x", Form("(%d -> %d) clust position", FOOT_N, POS), N_STRIPS, 0, N_STRIPS);
}

void TFOOTCalCont::Clean(Option_t* option) noexcept {
	this->N = 0; 
}

void TFOOTCalCont::AddCluster(double mean_x, double e, double mult, ClusterType ct) {
	if(N == static_cast<int>(CAPACITY)) {
		WARN("\'%s\' Trying to add a cluster, but capacity %zu reached. Rejecting it.\n", GetName(), CAPACITY);
		return;
	}
	_x[N] = mean_x;
	_e[N] = e;
	_m[N] = mult;
	_t[N] = ct;
	++N;
}

IMPL_CONTAINER_METHODS(TFOOTCalCont)

ClassImp(TFOOTCalCont);
