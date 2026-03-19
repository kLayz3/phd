#pragma once

#include "TFOOTMapCont.h"
#include "monad/monad.hxx"
#include "util/json_struct_def.hh"
#include "TFOOTCalCont.h"

class TH2I;

struct FOOTHit {
	double e,x,y,z;
	
	FOOTHit() = default;
	FOOTHit(double e, double x, double y, double z) :
		e(e), x(x), y(y), z(z) {}

	virtual ~FOOTHit() = default;
	ClassDef(FOOTHit, 1);
};


struct RNFOOTHit {
	std::vector<FOOTHit> hits;
	RNFOOTHit() = default;

	inline void Clean() noexcept { hits.clear(); }
	virtual ~RNFOOTHit() = default;
	ClassDef(RNFOOTHit, 1);
};

struct TFOOTHitCont : TContainer<RNFOOTHit> {
	inline static nlohmann::json setup {};
	inline static FOOTBoxParam _box;
	static constexpr int N_PAIRS = N_FOOT_DETECTORS / 2;
	
	FOOTBoxParam* box;
	std::string* setupName;
	TH2I* h_corr_all[N_PAIRS];
	TH2I* h_corr_all_sorted[N_PAIRS];
	TH2I* h_corr_gud[N_PAIRS];
	TH1I* h_single_all[N_FOOT_DETECTORS];
	TH1I* h_single_gud[N_FOOT_DETECTORS];
	TFOOTHitCont();

	void Setup() override;
	void Init(TDictInfo ) override;
};
