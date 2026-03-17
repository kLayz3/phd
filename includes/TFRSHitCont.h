#pragma once

#include "monad/monad.hxx"
#include "util/json_struct_def.hh"
#include <utility>

#include "Rtypes.h"
#include "TH2D.h"

struct FRSTargetParam {
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(double, width,     0.0, 0);
	ADD_SERIALIZABLE_FIELD(double, thickness, 0.0, 1);

	FRSTargetParam() = default;
	virtual ~FRSTargetParam() = default;
	ClassDef(FRSTargetParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSTargetParam, 1)

struct FRSIdParam {
	GET_HELP_AUX_IMPL;
	 
	ADD_SERIALIZABLE_FIELD(double,         dist,          0.0,  0); /* Distance relative to iron yoke at S2 (last quad b4 air). */
	ADD_SERIALIZABLE_FIELD(double,         magnification, 1.0,  1); 
	ADD_SERIALIZABLE_FIELD(double,         brho,          10.0, 2); /* Brho at the entrance to that focal point. */
	ADD_SERIALIZABLE_FIELD(double,         dispersion,    0.0,  3);
	ADD_SERIALIZABLE_FIELD(double,         zfocus,        1000, 4);
	ADD_SERIALIZABLE_FIELD(FRSTargetParam, target,        {},   5);

	FRSIdParam() = default;
	virtual ~FRSIdParam() = default;
	ClassDef(FRSIdParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSIdParam, 5)

struct RNFRSHit {
	/* At this level, FRS ID in a single focal point is determined uniquely by:
	 * - β velocity
	 * - (x,y) at target
	 * - (a,b) at target
	 * We **do not** allow multiple candidates for the incoming particle.
	 * It must either resolve to a single unique track or nothing. Multihit is anyway
	 * invalidated, as FOOT's can reliably map-out only a single reaction. */

	/* One per focal point. S2/S3 */
	struct Id {
		double x = NAN;
		double y = NAN;
		double a = NAN;
		double b = NAN;
		double beta = NAN;

		void Clean() noexcept {
			x=y=a=b=beta = NAN;
		}
		virtual ~Id() = default;
		ClassDef(Id, 1);
	} s2_bt, s2_at, s3; /* `bt` == before target; `at` == after target */
	
	inline void Clean() noexcept {
		s2_bt.Clean();
		s2_at.Clean();
		s3.Clean();
	}
	virtual ~RNFRSHit() = default;
	ClassDef(RNFRSHit, 1);
};

struct TFRSHitCont : TContainer<RNFRSHit> {
	inline static nlohmann::json setup {};
	inline static FRSIdParam _s2p, _s3p;
	inline static FRSTargetParam _sTar;

	TH2D *h2_track_x, *h2_track_y;

	FRSIdParam *s2p, *s3p;
	FRSTargetParam *sTar;
	std::string *setupName;
		
	TFRSHitCont();

	void Init(TDictInfo info) override;
	void Setup() override;
};
