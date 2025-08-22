#include "TFOOTPedestalCont.h"
#include "TContainer.h"
#include "TH2.h"
#include <algorithm>
#include <cmath>
#include <string>

TFOOTPedestalCont::TFOOTPedestalCont() {}
TFOOTPedestalCont::TFOOTPedestalCont(int N) : TContainer(Form("FOOT%d", N)), FOOT_N(N) {}
TFOOTPedestalCont::~TFOOTPedestalCont() {}

void TFOOTPedestalCont::Init(TDictInfo info) {
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
	
	assert(GetOwnedTOnceObjects().size() == 0 && "Don't call `SetId` twice. Clear the owned objects first.");

	h2_raw     = RegisterObject<TH2I>("h2_raw" , Form("Raw FOOT%s", GetName()), 640,0,640,4096,0,4096);
	h2_mid     = RegisterObject<TH2D>("h2_mid", Form("Raw FOOT%s - Global Pedestal", GetName()), 640,0,640,3000,-500,2500);
	h2_corr    = RegisterObject<TH2D>("h2_corr", Form("Corrected FOOT%s", GetName()), 640,0,640,6000,-500,2500);
	gr_s0      = RegisterObject<TGraph>("sigma0", N_ASIC);
	gr_s1      = RegisterObject<TGraph>("sigma1", N_ASIC);
	gped       = RegisterObject<std::array<double, N_STRIPS>>("gped", {});
	gped_s     = RegisterObject<std::array<double, N_STRIPS>>("gped_sigma", {});
	gped_sf    = RegisterObject<std::array<double, N_STRIPS>>("ped_sigma_corr", {});
	bad_strips = RegisterObject<std::vector<int>>("bad_strips", {});	
	h2_ped_off_med  = RegisterObject<TH2D>("h2_ped_off_med" , Form("Fine pedestal %s offset calculated via median", GetName()), 10,0,10, 1000,-100,100); 
	h2_ped_off_avg  = RegisterObject<TH2D>("h2_ped_off_avg" , Form("Fine pedestal %s offset calculated via trimmed average", GetName()), 10,0,10, 1000,-100,100);
	h2_ped_off_diff = RegisterObject<TH2D>("h2_ped_off_diff", Form("Fine pedestal %s offset diff (median - trimmed average)", GetName()), 10,0,10, 100,-10,10);

	gr_s0->SetMarkerStyle(22);
	gr_s0->SetMarkerSize(1);
	gr_s0->SetMarkerColor(kRed);
	gr_s1->SetMarkerStyle(23);
	gr_s1->SetMarkerSize(1);
	gr_s1->SetMarkerColor(kGreen);
}

void TFOOTPedestalCont::Clean(Option_t* option) noexcept {
	(void)option;
	std::fill_n(FOOTE, LEN(FOOTE), std::nan(""));
}

IMPL_CONTAINER_METHODS(TFOOTPedestalCont)

ClassImp(TFOOTPedestalCont);
