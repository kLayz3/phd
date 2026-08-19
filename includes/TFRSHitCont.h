#pragma once

#include "monad/monad.hxx"

#include "TFRSCalCont.h"
#include "util/json_struct_def.hh"
#include <utility>

#include "Rtypes.h"
#include "TH2D.h"

/* In this low-Z experiments, we always perform the PiD with
 * different scintillators, and not directly with gas detectors. */
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

struct FRSToFSingle {
    GET_HELP_AUX_IMPL;

    template<typename T>
    using A = std::array<T, 2>;
    
    ADD_SERIALIZABLE_FIELD(A<u32>,    combo, {}, 0);
    ADD_SERIALIZABLE_FIELD(A<double>, par,   {}, 1);

	FRSToFSingle() = default;
	virtual ~FRSToFSingle() = default;
	ClassDef(FRSToFSingle, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSToFSingle, 1)

struct FRSToFParam {
    GET_HELP_AUX_IMPL;

    using FRSToFSingleVec = std::vector<FRSToFSingle>;
    ADD_SERIALIZABLE_FIELD(FRSToFSingleVec, ToF, {}, 0);
	
    FRSToFParam() = default;
	virtual ~FRSToFParam() = default;
	ClassDef(FRSToFParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSToFParam, 0)

struct RNFRSHit {
	static constexpr i32 S2_BT_TRACKING_INCLUDE_SCI21_MASK            = 0x01;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC21_MASK            = 0x02;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC22_MASK            = 0x04;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC23_MASK            = 0x08;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_POINTLIKE_TARGET_MASK = 0x10;

	/* At this level, FRS ID in a single focal point is determined uniquely by:
     * - (A,Q) particle ID
	 * - β velocity
	 * - (x,y) at z=0, nominal
	 * - (a,b) angles
	 * We **do not** allow multiple candidates for the incoming particle.
	 * It must either resolve to a single unique track or nothing. Multihit is anyway
	 * invalidated, as FOOT's can reliably map-out only a single reaction. */

	/* S2/S3 */
	struct Id {
        f64 A  = NAN;
        f64 Q  = NAN;
		f64 x0 = NAN; // at nominal z = 0 ; FRS standard coordinates
		f64 y0 = NAN; // at nominal z = 0 ; FRS standard coordinates
		f64 ax = NAN;
		f64 ay = NAN;
		f64 beta = NAN;

		i64 code = 0; // some metadata fed from the processor.
		
		Id() = default;

		inline void Clean() noexcept {
			A=Q=x0=y0=ax=ay=beta = NAN;
			code = 0;
		}
		virtual ~Id() = default;
		ClassDef(Id, 1);
	};
	Id s2_bt, s2_at, s3; /* `bt` == before target; `at` == after target */
	
	double xT = NAN; // at EXPERT target
	double yT = NAN; // at EXPERT target

	std::string DecodeS2() const noexcept;
	std::string DecodeS3() const noexcept;
	
	/* To keep a bit of a debug-handle on previous step,
	 * just send the RNFRSCal here. */
	RNFRSCal cal;

	RNFRSHit() = default;
	inline void Clean() noexcept {
		s2_bt.Clean();
		s2_at.Clean();
		s3.Clean();
		xT = NAN; yT = NAN;
		//cal.Clean(); // dont need since it's just copied over from prev. step
	}

	virtual ~RNFRSHit() = default;
	ClassDef(RNFRSHit, 1);
};

struct TFRSHitCont : TContainer<RNFRSHit> {
	TH2D *h2_track_x, *h2_track_y;
	TH2D *h2_target_xy;

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param{};
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_param{};

	FRSIdParam *s2p, *s3p;
	FRSTargetParam *sTar;
    FRSToFParam *sTof;
	std::string *setupName;
	
	TFRSHitCont();

	void Init(TDictInfo info) override;
	void Setup() override;
};

namespace mnd::geom { struct Line3D; }
mnd::geom::Line3D RNTrackToLine3D(const RNFRSHit::Id& );
