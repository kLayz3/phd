#include "TFOOTCalCont.h"
#include "TH1I.h"
#include "json_struct_def.hh"
#include "nlohmann/json.hpp"
using json = nlohmann::json;

RNFOOTCluster::RNFOOTCluster(double x, double e, double m, ClusterType t) :
	fCX(x), fCE(e), fCM(m), fCT(t) {}

/* ------------------------------------------------------- */

RNFOOTCal::RNFOOTCal() { 
	fCl.reserve(INIT_CAPACITY); 
	_fBadE.reserve(N_STRIPS);
	_fHeClSize1.reserve(N_STRIPS);
}
void RNFOOTCal::Clean() noexcept {
	fCl.clear(); 
	_fBadE.clear();
	_fHeClSize1.clear();
}
std::vector<double> RNFOOTCal::E() const noexcept {
	std::vector<double> res;
	res.reserve(fCl.size());
	for(auto const& c : fCl) res.push_back(c.fCE);
	return res;
}
std::vector<double> RNFOOTCal::X() const noexcept {
	std::vector<double> res;
	res.reserve(fCl.size());
	for(auto const& c : fCl) res.push_back(c.fCX);
	return res;
}

/* ------------------------------------------------------- */

/* Initialized with 2 keys: "FOOT_ID" and "Setup" */
void TFOOTCalCont::Init(TDictInfo info) {
	constexpr const char* id_key = "ID";
	auto it = info.find(id_key);
	if(it == info.end())
		ERROR("\"%s\" key not found in the info hashmap.", id_key);
	try {
		FOOT_N = std::stoi(it->second);
	} catch(const std::exception& e) {
		ERROR("For key \"%s\", value found: " EMPH(%s) " but unparsable to integer. Err: %s", id_key, it->second.c_str(), e.what());
	}

	it = info.find("Setup");	
	if(it == info.end())
		ERROR("TFOOTCalCont::Init(): \"Setup\" key not found for info (%s).\n" 
			, mnd::type_name<TDictInfo>().c_str());
	json j {};
	const std::string& file_name = it->second;
	auto f = mnd::get_maybe_ifstream(file_name);
	if(!f.has_value())
		ERROR("File \'%s\' not found or not openable.\n", file_name.c_str());
	j = json::parse(std::move( f.value() ));
	auto j_it = j.find(Form("FOOT%d", FOOT_N));
	if(j_it == j.end()) 
		ERROR("Trying to set up FOOT[%d] cal, but setup file doesn't contain the key \"FOOT%d\".", FOOT_N, FOOT_N);
	
	json& jf = j_it.value();
	UNROLL_JSON_PARAM(par, jf, 7);

	this->SetName(Form("FOOT%d", par.N));
}

using A2 = std::array<double, 2>;
template<> void Add(A2& lhs, const A2& rhs) {}

void TFOOTCalCont::Setup() {
	if(strlen(GetName()) == 0) ERROR("Setup called before Init?");
	h1_raw_mult = RegisterObject<TH1I>("h1_raw_mult", Form("FOOT(%2d: %d) raw multiplicity", FOOT_N, par.N), 50,0,50);
	h1_mult     = RegisterObject<TH1I>("h1_mult", Form("FOOT (%2d: %d) multiplicity", FOOT_N, par.N), 200,0,50);
	h1_dE       = RegisterObject<TH1I>("h1_dE", Form("FOOT(%2d: %d) energy", FOOT_N, par.N), 5000,0,1000);
	h1_X        = RegisterObject<TH1I>("h1_x", Form("FOOT(%2d: %d) clust position", FOOT_N, par.N), 5*N_STRIPS, 0, N_STRIPS);
	h1_cl_type  = RegisterObject<TH1I>("h1_cl_type", Form("(%2d:%d) clust type", FOOT_N, par.N), 5, 0, 5);
	h1_dE_m1    = RegisterObject<TH1I>("h1_dE_m1", Form("(%2d:%d) energy with multiplicity 1", FOOT_N, par.N), 2500, 0, 500);
	h1_dE_m2    = RegisterObject<TH1I>("h1_dE_m2", Form("(%2d:%d) energy with multiplicity 2", FOOT_N, par.N), 2500, 0, 500);
	h1_dE_m3    = RegisterObject<TH1I>("h1_dE_m3", Form("(%2d:%d) energy with multiplicity 3", FOOT_N, par.N), 2500, 0, 500);
	h1_sn_ratio = RegisterObject<TH1I>("h1_sn_ratio", Form("(%2d:%d) ratio neighbouring vs. seed value (mult <= 3)", FOOT_N, par.N), 500, 0, 5);
	
	setup = RegisterObject<FOOTParam>("setup", mnd::noop_fn<FOOTParam>(), this->par /* copy ctor */);
}

ClassImp(RNFOOTCluster);
ClassImp(RNFOOTCal);

ClassImp(FOOTParam);
