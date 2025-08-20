#include "TFOOTPedestalCont.h"
#include "AuxFunctions.hh"
#include "TContainer.h"
#include "TH2.h"
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
	
	assert(sizeof(GetOwnedTOnceObjects()) == 0 && "Don't call `SetId` twice. Clear the owned objects first.");

	h2_raw  = RegisterObject<TH2I>(Form("%s_h2_raw", GetName()), Form("Raw FOOT%s", GetName()), 640,0,640,4096,0,4096);
	h2_corr = RegisterObject<TH2D>(Form("%s_h2_corr", GetName()), Form("Corrected FOOT%s", GetName()), 640,0,640,3000,-500,2500);
	h2_s0   = RegisterObject<TH2D>(Form("%s_h2_sigma0", GetName()), Form("Sigma0 FOOT%s", GetName()), 640,0,640,100,0,20);
	h2_s1   = RegisterObject<TH2D>(Form("%s_h2_sigma1", GetName()), Form("Sigma1 FOOT%s", GetName()), 640,0,640,100,0,20);
	gped    = RegisterObject<std::array<double, N_STRIPS>>(Form("%s_gped", GetName()), {});
	gped_s  = RegisterObject<std::array<double, N_STRIPS>>(Form("%s_gped_s", GetName()), {});

	h2_ped_off_med = RegisterObject<TH2D>(Form("%s_h2_ped_off_med", GetName()), Form("Fine pedestal offset calculated via median"), 10,0,10, 1000,-100,100); 
	h2_ped_off_avg = RegisterObject<TH2D>(Form("%s_h2_ped_off_avg", GetName()), Form("Fine pedestal offset calculated via trimmed average"), 10,0,10, 1000,-100,100);
	h2_ped_off_diff = RegisterObject<TH2D>(Form("%s_h2_ped_off_diff", GetName()), Form("Fine pedestal offset diff (median - trimmed average)"), 10,0,10, 100,-10,10);
}

void TFOOTPedestalCont::Clean(Option_t* option) noexcept {
	(void)option;
	has_data = false;
}

IMPL_CONTAINER_SETUP(TFOOTPedestalCont)
ClassImp(TFOOTPedestalCont);
