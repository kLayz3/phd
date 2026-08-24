#include "TFOOTCalCont.h"
#include "util/JSONParser.h"
#include "nlohmann/json.hpp"
#include "util/json_struct_def.hh"
#include "util/FFT.h"
#include "filesystem"

using json = nlohmann::json;
namespace fs = std::filesystem;

double FOOTDeltaFFT::Evaluate(const double x) const {
	return FFTW::Evaluate(this->c, this->n, x, -0.5, 0.5); 	
}

double FOOTParam::E(const RNFOOTCluster& clust) const noexcept {
	double ce = clust.fCE;
	double d  = clust.Delta();
	double cx = clust.fCX;

	ce *= this->gain.CorrectionFactor(cx, ce);
	ce /= this->de.CorrectionFactor(d);

	return ce;
};

double FOOTParam::Q(const RNFOOTCluster& clust) const noexcept {
	double e = this->E(clust); 
	return this->Q(e);
}

double FOOTParam::BarePosition(const RNFOOTCluster& clust) const noexcept {
	return R() * (clust.fCX - DETECTOR_MIDPOINT) * STRIP_TO_MM;
}
double FOOTParam::X0(const RNFOOTCluster& clust) const noexcept {
	return BarePosition(clust) + delta_p;
}

RNFOOTCluster::RNFOOTCluster(double x, double e, u32 m, ClusterType t, FOOTClusterFit fit) :
	fCX(x), fCE(e), fCM(m), fCT(t), fit{fit} {}

/* ------------------------------------------------------- */

RNFOOTCal::RNFOOTCal() {
	fCl.reserve(INIT_CAPACITY);
	fRaw.reserve(N_STRIPS);
}
void RNFOOTCal::Clean() noexcept {
	fCl.clear();
	fRaw.clear();
    w_ev_type = EventType::Undetermined;
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

static std::string full_setup_file_name {};

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
		
		UNROLL_JSON_PARAM(bpar, j["box"], 7);
		
		/* `file_name` could've been passed relative. Fetch the full path and bake it in. */
		try {
			fs::path pfull = fs::canonical(file_name);
			full_setup_file_name = pfull.c_str();
		} catch(std::filesystem::filesystem_error const& e) {
			ERROR("Fetching full path from: \'%s\' failed? Reason: %s", file_name.c_str(), e.what());
		}
	}
}

void TFOOTCalCont::Setup() {
	if(strlen(GetName()) == 0) ERROR("Setup called before Init() ? Or name field set to '' ?");
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

	setup = RegisterObject<FOOTParam>("setup", this->par /* copy ctor */);
	
	/* In the cal step, box object not used. It's just written down for helper scripts to process stuff. */
	if(should_register_box_) {
		box = RegisterObject<FOOTBoxParam>("box", this->bpar);
		setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), full_setup_file_name); 
	}
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
ClassImp(ExpertTarget);
ClassImp(FOOTBoxParam);
