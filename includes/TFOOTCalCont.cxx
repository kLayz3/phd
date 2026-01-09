#include "TFOOTCalCont.h"
#include "TH1I.h"
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
	
	auto it = info.find("Setup");	
	if(it == info.end()) {
		WARN("TFOOTCalCont::Init(): `Setup` key not found for info (%s). " 
			"Is fine, will default to (%.2f, %.2f) thresholds (in units of sigma).\n", mnd::type_name<TDictInfo>().c_str(),
			TFOOTCalCont::CENTRE_THR_DEFAULT, TFOOTCalCont::NEIGHB_THR_DEFAULT);
	}
	else {
		json j {};
		const std::string& file_name = it->second;
		auto f = mnd::get_maybe_ifstream(file_name);
		if(!f.has_value())
			ERROR("File \'%s\' not found or not openable.\n", file_name.c_str());

		try {
			j = json::parse(std::move( f.value() ));
			auto j_it = j.find(Form("FOOT%d", FOOT_N));
			if(j_it == j.end()) 
				ERROR("FOOT[%d -> %d], setup file doesn't contain the key \"FOOT%d\".", FOOT_N, POS, FOOT_N);
			
			json& jf = j_it.value();
			if(j_it = jf.find("c_threshold"); j_it != jf.end()) {
				c_threshold = j_it.value().get<double>();
			}
			if(j_it = jf.find("n_threshold"); j_it != jf.end()) {
				n_threshold = j_it.value().get<double>();
			}
			assert(c_threshold > n_threshold && Form("FOOT[%d -> %d] cthreshold(%.2f) < nthreshold(%.2f). Not allowed,",
				FOOT_N, POS, c_threshold, n_threshold));

		} catch(std::exception const& e) {
			ERROR("Setup parse failed. Reason: %s\n", e.what());
		}
	}

	this->SetName(Form("FOOT%d", POS));
}

using A2 = std::array<double, 2>;
template<> void Add(A2& lhs, const A2& rhs) {}

void TFOOTCalCont::Setup() {
	if(strlen(GetName()) == 0) ERROR("Setup called before Init?");
	h1_raw_mult = RegisterObject<TH1I>("h1_raw_mult", Form("FOOT(%2d: %d) raw multiplicity", FOOT_N, POS), 50,0,50);
	h1_mult     = RegisterObject<TH1I>("h1_mult", Form("FOOT (%2d: %d) multiplicity", FOOT_N, POS), 200,0,50);
	h1_dE       = RegisterObject<TH1I>("h1_dE", Form("FOOT(%2d: %d) energy", FOOT_N, POS), 5000,0,1000);
	h1_X        = RegisterObject<TH1I>("h1_x", Form("FOOT(%2d: %d) clust position", FOOT_N, POS), 5*N_STRIPS, 0, N_STRIPS);
	h1_cl_type  = RegisterObject<TH1I>("h1_cl_type", Form("(%2d:%d) clust type", FOOT_N, POS), 5, 0, 5);
	h1_dE_m1    = RegisterObject<TH1I>("h1_dE_m1", Form("(%2d:%d) energy with multiplicity 1", FOOT_N, POS), 2500, 0, 500);
	h1_dE_m2    = RegisterObject<TH1I>("h1_dE_m2", Form("(%2d:%d) energy with multiplicity 2", FOOT_N, POS), 2500, 0, 500);
	h1_dE_m3    = RegisterObject<TH1I>("h1_dE_m3", Form("(%2d:%d) energy with multiplicity 3", FOOT_N, POS), 2500, 0, 500);
	h1_sn_ratio = RegisterObject<TH1I>("h1_sn_ratio", Form("(%2d:%d) ratio neighbouring vs. seed value (mult <= 3)", FOOT_N, POS), 500, 0, 5);
	
	threshold = RegisterObject<A2>("threshold", mnd::noop_fn<A2>(), {c_threshold, n_threshold});
}

ClassImp(RNFOOTCluster);
ClassImp(RNFOOTCal);
