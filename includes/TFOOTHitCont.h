#pragma once

#include "monad/monad.hxx"
#include "json_struct_def.hh"
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

struct FOOTBoxParam {
	GET_HELP_AUX_IMPL
	using T1 = 	std::array<double, 4>;
	ADD_SERIALIZABLE_FIELD(double, z0,          NAN, 0);
	ADD_SERIALIZABLE_FIELD(double, width_inner, NAN, 1);
	ADD_SERIALIZABLE_FIELD(double, width_outer, NAN, 2);
	ADD_SERIALIZABLE_FIELD(double, dx,          NAN, 3);
	ADD_SERIALIZABLE_FIELD(double, dy,          NAN, 4);
	ADD_SERIALIZABLE_FIELD(double, da,          NAN, 5);
	ADD_SERIALIZABLE_FIELD(double, db,          NAN, 6);
	virtual ~FOOTBoxParam() = default;
	ClassDef(FOOTBoxParam, 1);
};
ADD_STD_TYPE_RESOLUTION_(FOOTBoxParam, 6)

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
