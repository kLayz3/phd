#include "TFOOTMapCont.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TGraph.h"
#include <cmath>
#include <string>

TFOOTMapCont::TFOOTMapCont(int N) : TContainer(Form("FOOT%d", N)), FOOT_N(N) {}

void TFOOTMapCont::Init(TDictInfo info) {
	auto n_it = info.find("FOOT_ID");
	if(n_it == info.end())
		ERROR("FOOT_ID key not found in the info hashmap.");
	i32 n;
	try {
		n = std::stoi(n_it->second);
	} catch(const std::exception& e) {
		ERROR("FOOT_ID found: " EMPH(%s) " but unparsable to integer. Err: %s", n_it->second.c_str(), e.what());  
	}
	FOOT_N = n;
	this->SetName(Form("FOOT%d", n));
}

void TFOOTMapCont::Setup() {
	if(FOOT_N < 0) ERROR("FOOT index not set. It is %d. Did you call `Init` ?", FOOT_N);

	h2_raw_tmp = RegisterObject<TH2I>("h2_raw_tmp" , Form("Raw FOOT%d Temp", FOOT_N), 640,0,640,4096,0,4096);
	h2_raw     = RegisterObject<TH2I>("h2_raw" , Form("Raw FOOT%d", FOOT_N), 640,0,640,4096,0,4096);
	h2_corr    = RegisterObject<TH2D>("h2_corr", Form("Corrected FOOT%d", FOOT_N), 640,0,640,6000,-500,2500);

	h2_ped0  = RegisterObject<TH2D>("h2_ped0", Form("Initial pedestal FOOT%d per strip", FOOT_N), 
		N_STRIPS, 0, N_STRIPS, 200,0, 10);
	h2_sigma0  = RegisterObject<TH2D>("h2_sigma0", Form("Initial sigma FOOT%d per strip", FOOT_N), 
		N_STRIPS, 0, N_STRIPS, 200,0, 10);
	gr_s1      = RegisterObject<TGraph>("sigma1", N_ASIC);
	gped       = RegisterObject<std::array<double, N_STRIPS>>("gped", {});
	gped_s     = RegisterObject<std::array<double, N_STRIPS>>("gped_sigma", {});
	gped_sf    = RegisterObject<std::array<double, N_STRIPS>>("ped_sigma_corr", {});
	bad_strips = RegisterObject<std::vector<int>>("bad_strips", {});	
	h2_ped_off_med  = RegisterObject<TH2D>("h2_ped_off_med" , Form("Fine pedestal FOOT%d offset calculated via median", FOOT_N), 10,0,10, 1000,-100,100); 
	h2_ped_off_avg  = RegisterObject<TH2D>("h2_ped_off_avg" , Form("Fine pedestal FOOT%d offset calculated via trimmed average", FOOT_N), 10,0,10, 1000,-100,100);
	h2_ped_off_diff = RegisterObject<TH2D>("h2_ped_off_diff", Form("Fine pedestal FOOT%d offset diff (median - trimmed average)", FOOT_N), 10,0,10, 100,-10,10);

	gr_s1->SetMarkerStyle(23);
	gr_s1->SetMarkerSize(1);
	gr_s1->SetMarkerColor(kGreen);
}

ClassImp(RNFOOTMap);
