#include "TFOOTCalCont.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TParameter.h"
#include "util/JSONParser.h"
#include "nlohmann/json.hpp"
#include "util/json_struct_def.hh"
#include "util/FFT.h"
using json = nlohmann::json;

double FOOTDeltaFFT::Evaluate(const double x) const {
	return FFTW::Evaluate(this->c, this->n, x, -0.5, 0.5); 	
}

RNFOOTCluster::RNFOOTCluster(double x, double e, u32 m, ClusterType t, FOOTClusterFit fit) :
	fCX(x), fCE(e), fCM(m), fCT(t), fit{fit} {}

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

/* Initialized with 2 keys: "ID" and "Setup" */
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
	const std::string& file_name = it->second;
	json j = ParseJSON(file_name);
	
	auto j_it = j.find(Form("FOOT%d", FOOT_N));
	if(j_it == j.end()) 
		ERROR("Trying to set up FOOT[%d] cal, but setup file doesn't contain the key \"FOOT%d\".", FOOT_N, FOOT_N);
	
	json& jf = j_it.value();
	UNROLL_JSON_PARAM(par, jf, 9);

	par.de10_index_ = FOOT_N;
	if(par.N == -1) 
		ERROR("Parsed the setup file fine, but the \"N\" table entry for FOOT%d not found. "
			"It is mandatory to label the FOOT's!", FOOT_N);
	/* Note, if two FOOT's in the setup have identical `.N` field, it will throw a
	 * wildest RNTuple error 'cause of columns of identical name... */
	this->SetName(Form("FOOT%d", par.N));

	if(should_register_box_) {
		if(!j.contains("box"))
			ERROR("Container meant to also register the FOOT box, but inside %s file, \"box\" key not found?",
				file_name.c_str());
		
		UNROLL_JSON_PARAM(bpar, j["box"], 8);
	}
}

using A2 = std::array<double, 2>;
template<> void Add(A2& lhs, const A2& rhs) {}

void TFOOTCalCont::Setup() {
	if(strlen(GetName()) == 0) ERROR("Setup called before Init? or name field set to '' ?");
	h1_mult     = RegisterObject<TH1I>("h1_mult", Form("FOOT (%2d: %d) multiplicity", FOOT_N, par.N), 200,0,50);
	h1_dE       = RegisterObject<TH1I>("h1_dE", Form("FOOT(%2d: %d) energy", FOOT_N, par.N), 5000,0,1000);
	h1_X        = RegisterObject<TH1I>("h1_x", Form("FOOT(%2d: %d) clust position", FOOT_N, par.N), 5*N_STRIPS, 0, N_STRIPS);
	h1_cl_type  = RegisterObject<TH1I>("h1_cl_type", Form("(%2d:%d) clust type", FOOT_N, par.N), 5, 0, 5);
	h1_dE_m1    = RegisterObject<TH1I>("h1_dE_m1", Form("(%2d:%d) energy with multiplicity 1", FOOT_N, par.N), 2500, 0, 500);
	h1_dE_m2    = RegisterObject<TH1I>("h1_dE_m2", Form("(%2d:%d) energy with multiplicity 2", FOOT_N, par.N), 2500, 0, 500);
	h1_dE_m3    = RegisterObject<TH1I>("h1_dE_m3", Form("(%2d:%d) energy with multiplicity 3", FOOT_N, par.N), 2500, 0, 500);
	h1_cl_sigma = RegisterObject<TH1I>("h1_cl_sigma", Form("(%2d:%d) fitted cluster width", FOOT_N, par.N), 500, 0, 5);
	h2_mult_e = RegisterObject<TH2I>("h2_mult_e", Form("FOOT(%2d:%d) multiplicity and energy", FOOT_N, par.N),
			2000,0,2000, 10, 1, 10);
	h2_mult_e->GetXaxis()->SetTitle("Individual energy per strip hit");
	h2_mult_e->GetYaxis()->SetTitle("Multiplicity");

	setup = RegisterObject<FOOTParam>("setup", mnd::noop_fn<FOOTParam>(), this->par /* copy ctor */);
	
	/* In the cal step, box object not used. It's just written down for helper scripts to process stuff. */
	if(should_register_box_)
		box = RegisterObject<FOOTBoxParam>("box", mnd::noop_fn<FOOTBoxParam>(), this->bpar);
}

ClassImp(FOOTClusterFit);
ClassImp(RNFOOTCluster);
ClassImp(RNFOOTCal);
ClassImp(FMultiPoly);
ClassImp(FOOTAsicGainParam);
ClassImp(FOOTReferentADCMeasurement);
ClassImp(FOOTGainParam);
ClassImp(FOOTDeltaFFT);
ClassImp(FOOTDeltaParam);
ClassImp(FOOTParam);
ClassImp(FOOTBoxParam);
