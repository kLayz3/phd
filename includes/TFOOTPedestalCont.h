#pragma once

#include "TContainer.hxx"
#include "libs.hh"
#include <array>
#include <vector>
#include "TOnce.hxx"

#define _FOOT_N_STRIPS         640
#define _FOOT_N_ASIC            10
#define _FOOT_N_STRIPS_PER_ASIC 64

class TH2I;
class TH2D;
class TGraph;

struct RNFOOTPedestalCont {	
	double FOOTE[_FOOT_N_STRIPS];
	inline void Clean() noexcept { std::fill_n(FOOTE, _FOOT_N_STRIPS, std::nan("")); }
	ClassDef(RNFOOTPedestalCont, 1);
};

class TFOOTPedestalCont : public TContainer<RNFOOTPedestalCont> {
public:
	static constexpr int N_STRIPS = _FOOT_N_STRIPS;
	static constexpr int N_ASIC = _FOOT_N_ASIC;
	static constexpr int N_STRIPS_PER_ASIC = _FOOT_N_STRIPS_PER_ASIC;

public: // Inputs from Go4.
	u32* _FOOT;
	u32* _FOOTE; // is array: [N_STRIPS]
	i32 FOOT_N;

public:
	// Objects held by the TContainer::_vc  
	TH2I* h2_raw; 
	TH2D* h2_mid; 
	TH2D* h2_corr;
	TGraph* gr_s0;
	TGraph* gr_s1;
	std::array<double, N_STRIPS> *gped, *gped_s, *gped_sf; // indices [0], [1], [2], ..., [639]
	std::vector<int> *bad_strips; // [0, 1, 2, ..., 639]
	TH2D* h2_ped_off_med;
	TH2D* h2_ped_off_avg;
	TH2D* h2_ped_off_diff;

public:
	TFOOTPedestalCont(int);
	
	void Init(TDictInfo info) /* override */;
};
