#include "TFOOTHitCont.h"
#include "TH2I.h"
#include "TH1I.h"
#include "util/JSONParser.h"

using json = nlohmann::json;

TFOOTHitCont::TFOOTHitCont() : TContainer("FOOT") {}

void TFOOTHitCont::Init(TDictInfo info) {
auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	
	setup = ParseJSON(file_name);
	setup["file_name"] = file_name;
	
	if(!setup.contains("box")) ERROR("\'box\' key not found in JSON file: %s\n", file_name.c_str());
		
	UNROLL_JSON_PARAM(_box, setup["box"], 7);
}

void Add(FOOTBoxParam&, const FOOTBoxParam&) {}
void TFOOTHitCont::Setup() {
	box = RegisterObject<FOOTBoxParam>("Box");
	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());
	for(int i=0; i<N_PAIRS; ++i) {
		h_corr_all[i] = RegisterObject<TH2I>(
			Form("h2_corr_all%d", i), 
			Form("Correlation in dE FOOT%d:FOOT%d for all clusters", 2*i, 2*i+1),
			300,0,3000,300,0,3000
		);
		h_corr_gud[i] = RegisterObject<TH2I>(
			Form("h2_corr_gud%d", i), 
			Form("Correlation in dE FOOT%d:FOOT%d for only paired up clusters", 2*i, 2*i+1),
			300,0,3000,300,0,3000
		);
		h_corr_all_sorted[i] = RegisterObject<TH2I>(
			Form("h2_corr_all_sorted%d", i), 
			Form("Correlation in dE FOOT%d:FOOT%d for sorted clusters (by energy)", 2*i, 2*i+1),
			300,0,3000,300,0,3000
		);
	}

	for(int i=0; i<N_FOOT_DETECTORS; ++i) {
		h_single_all[i] = RegisterObject<TH1I>(
			Form("h1_single_all%d", i), 
			Form("dE FOOT%d (after eta-correction) for all clusters", i),
			300,0,3000
		);
		h_single_gud[i] = RegisterObject<TH1I>(
			Form("h1_single_gud%d", i), 
			Form("dE FOOT%d (after eta-correction) for only paired-up clusters", i),
			300,0,3000
		);
	}
}

ClassImp(FOOTHit);
ClassImp(RNFOOTHit);
