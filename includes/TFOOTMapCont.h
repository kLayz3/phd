#pragma once

#include "monad/monad.hxx"
#include "Scaler.h"

#include <array>
#include <vector>

#define N_FOOT_DETECTORS 8

#define _FOOT_N_STRIPS         640
#define _FOOT_N_ASIC            10
#define _FOOT_N_STRIPS_PER_ASIC 64

class TH2I;
class TH2D;
class TGraph;
template<typename T> class TParameter;

struct RNFOOTMap {
	static constexpr int N_STRIPS          = _FOOT_N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = _FOOT_N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = _FOOT_N_STRIPS_PER_ASIC; /* 64  */

	Scaler<32> timing; // nullable according to `RNFOOTMap::HasData()` predicate 
	std::array<double, N_STRIPS> FOOTE;
	std::array<double, N_ASIC> common_offset;

	inline void Clean() noexcept { 
		FOOTE.fill( std::nan("") ); 
		common_offset.fill( std::nan("") ); 
	}

	inline bool HasData() const noexcept { return std::isfinite(FOOTE[0]); }
	inline mnd::Maybe<u32> T() const noexcept {
		return (HasData() and timing.initialized_) ? mnd::Maybe<u32>{timing.curr_data} : mnd::None;
	}
	
	inline mnd::Maybe<u64> TotalT() const noexcept {
		return HasData() ? mnd::Maybe<u64>{timing.cumulative} : mnd::None;
	}

	inline mnd::Maybe<u32> DeltaT() const noexcept {
		return (HasData() and timing.initialized_) ? mnd::Maybe<u32>{timing.increment} : mnd::None;
	}

	RNFOOTMap() = default;
	virtual ~RNFOOTMap() = default;
	ClassDef(RNFOOTMap, 1);
};

struct TFOOTMapCont : TContainer<RNFOOTMap> {
	static constexpr int N_STRIPS          = RNFOOTMap::N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = RNFOOTMap::N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = RNFOOTMap::N_STRIPS_PER_ASIC; /* 64  */
	static constexpr int N_GPED_CHANGE_TOLERANCE  = 5; /* Relative to intial gped, consecutive calculations can go +- 5 ADC units. */
	static constexpr int N_GPED_BINS_IN_TOLERANCE = 20; /* How many bins does the tolerance window, +- 5 span. */ 

	i32 FOOT_N = -1;
	
	TH2I* h2_raw_tmp; // Temporary container.
	TH2I* h2_raw;  // Total raw (from Go4) ADC spectrum.
	TH2D* h2_corr; // Total calibrated ADC spectrum.
	
	TH2D* h2_gped;       // Global pedestal calculated per strip; filled per batch.
	TH2D* h2_gped_sigma; // Global pedestal width calculated per strip; filled per batch. 

	TH2D* h2_gped_per_batch; // How global pedestal of each strip changes per batch.

	std::array<double, N_STRIPS> *ped_s; // indices [0], [1], [2], ..., [639]
	std::vector<int> *bad_strips; // [0, 1, 2, ..., 639]
	
	TH2D* h2_ped_off_med;
	TH2D* h2_ped_off_avg;
	TH2D* h2_ped_off_diff;
	TH1I* h1_entry_dt;

	TParameter<double>* dt_veto;

	void Setup() override;
	void Init(TDictInfo ) override;

	TFOOTMapCont(int);
	TFOOTMapCont() = default;
};
