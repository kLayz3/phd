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
	i32 FOOT_N;
	TFOOTMapCont() = default;
	TFOOTMapCont(int );

	void Init(TDictInfo ); /* override. */
};
