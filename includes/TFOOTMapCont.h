#pragma once

#include "core/TContainer.hxx"
#include "core/libs.hh"
#include "core/TOnce.hxx"

#include <array>
#include <vector>

#define _FOOT_N_STRIPS         640
#define _FOOT_N_ASIC            10
#define _FOOT_N_STRIPS_PER_ASIC 64

class TH2I;
class TH2D;
class TGraph;

struct RNFOOTMap {	
	double FOOTE[_FOOT_N_STRIPS]{};
	inline void Clean() noexcept { std::fill_n(FOOTE, _FOOT_N_STRIPS, std::nan("")); }
	
	RNFOOTMap() = default;
	
	virtual ~RNFOOTMap() = default;
	ClassDef(RNFOOTMap, 1);
};

struct TFOOTMapCont : TContainer<RNFOOTMap> {
	static constexpr int N_STRIPS          = _FOOT_N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = _FOOT_N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = _FOOT_N_STRIPS_PER_ASIC; /* 64  */
	
	i32 FOOT_N = -1;
	
	TH2I* h2_raw_tmp;
	TH2I* h2_raw; 
	TH2D* h2_corr;
	
	TH2D* h2_ped0;
	TH2D* h2_sigma0;

	TGraph* gr_s1;
	std::array<double, N_STRIPS> *gped, *gped_s, *gped_sf; // indices [0], [1], [2], ..., [639]
	std::vector<int> *bad_strips; // [0, 1, 2, ..., 639]
	TH2D* h2_ped_off_med;
	TH2D* h2_ped_off_avg;
	TH2D* h2_ped_off_diff;

	void Setup() override;
	void Init(TDictInfo ) override;

	TFOOTMapCont(int);
	TFOOTMapCont() = default;
};
