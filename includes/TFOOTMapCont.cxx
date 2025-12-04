#include "TFOOTMapCont.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TGraph.h"
#include <cmath>

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

	h2_raw_tmp = RegisterObject<TH2I>("h2_raw_tmp" , Form("Raw FOOT%d Temp", FOOT_N), N_STRIPS, 0, N_STRIPS, 4096,0,4096);
	h2_raw     = RegisterObject<TH2I>("h2_raw" , Form("Raw FOOT%d", FOOT_N), N_STRIPS, 0, N_STRIPS, 4096,0,4096);
	h2_corr    = RegisterObject<TH2D>("h2_corr", Form("Corrected FOOT%d", FOOT_N), N_STRIPS, 0, N_STRIPS, 6000,-500,2500);

	h2_gped  = RegisterObject<TH2D>("h2_gped", Form("Initial pedestal FOOT%d per strip", FOOT_N), 
		N_STRIPS, 0, N_STRIPS, 1000,0,1000);
	h2_gped_sigma  = RegisterObject<TH2D>("h2_gped_sigma", Form("Initial sigma FOOT%d per strip", FOOT_N), 
		N_STRIPS, 0, N_STRIPS, 200,0,10);

	h2_gped_per_batch  = RegisterObject<TH2D>("h2_gped_per_batch", Form("Global pedestal FOOT%d per data chuck batch", FOOT_N), 
		50, 0, 50,
		N_GPED_BINS_IN_TOLERANCE * N_STRIPS, -N_GPED_CHANGE_TOLERANCE, 2*N_GPED_CHANGE_TOLERANCE*N_STRIPS - N_GPED_CHANGE_TOLERANCE);

	ped_s     = RegisterObject<std::array<double, N_STRIPS>>("ped_sigma", {});
	bad_strips = RegisterObject<std::vector<int>>("bad_strips", {});	

	h2_ped_off_med  = RegisterObject<TH2D>("h2_ped_off_med" , Form("Fine pedestal FOOT%d offset calculated via median", FOOT_N), 10,0,10, 1000,-100,100); 
	h2_ped_off_avg  = RegisterObject<TH2D>("h2_ped_off_avg" , Form("Fine pedestal FOOT%d offset calculated via trimmed average", FOOT_N), 10,0,10, 1000,-100,100);
	h2_ped_off_diff = RegisterObject<TH2D>("h2_ped_off_diff", Form("Fine pedestal FOOT%d offset diff (median - trimmed average)", FOOT_N), 10,0,10, 100,-10,10);

}

ClassImp(RNFOOTMap);
