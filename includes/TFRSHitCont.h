#pragma once

#include "monad/monad.hxx"

#include "TFRSCalCont.h"
#include "util/RunsheetParser.h"
#include "util/json_struct_def.hh"
#include <utility>

#include "Rtypes.h"
#include "TH2D.h"

namespace phy { struct Nucleus; }
namespace mnd { struct RunsheetState; }

struct FRSToFSingle {
    GET_HELP_AUX_IMPL;

    template<typename T>
    using A = std::array<T, 2>;

    ADD_SERIALIZABLE_FIELD(A<u32>,    combo, {}, 0);
    ADD_SERIALIZABLE_FIELD(A<double>, par,   {}, 1);

	/* Calculate β based on TDC timing difference between given {par[0], par[1]} SCI's */
	double Beta(double dt) const noexcept;

	/* Calculate β based on entire cal step. It will just dereference sci containers with
	 * correct indices {par[0], par[1]} and calculate the velocity based on that.
	 * Always takes first hit, and have that as the reference. */
	double Beta(const RNFRSCal& ) const noexcept;

	/* Checks if the parsed paremeters from JSON are valid. */
	bool IsValid() const noexcept;

	FRSToFSingle() = default;
	virtual ~FRSToFSingle() = default;
	ClassDef(FRSToFSingle, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSToFSingle, 1)

struct FRSToFParam {
    GET_HELP_AUX_IMPL;

    using FRSToFVec = std::vector<FRSToFSingle>;
    ADD_SERIALIZABLE_FIELD(FRSToFVec, ToF, {}, 0);

	FRSToFSingle const* Get(u32, u32) const noexcept;

    FRSToFParam() = default;
	virtual ~FRSToFParam() = default;
	ClassDef(FRSToFParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FRSToFParam, 0)

/* This parameter does not come from setup file, but rather from the runsheet row. */
struct BeamInfo {
    GET_HELP_AUX_IMPL;
	using BrhoSeq = std::array<double, 4>;
	ADD_SERIALIZABLE_FIELD(u32,    A0,   0, 0); // Nucleon number (selected) coming into S2.
	ADD_SERIALIZABLE_FIELD(u32,    Z0,   0, 1); // Atomic  number (selected) coming into S2.
	ADD_SERIALIZABLE_FIELD(u32,    Q0,   0, 2); // Charge  state (selected) coming into S2.
	
	ADD_SERIALIZABLE_FIELD(BrhoSeq, R,  {}, 3); // Magnetic rigidities

	/* Retrieve the mean Bρ value at the given scintillator. */
	inline constexpr double const& BrhoAt(u32 sci_index) const noexcept {
		switch(sci_index) {
			case RNFRSCal::SCI21_I: return R[1];
			case RNFRSCal::SCI22_I: return R[2];
			case RNFRSCal::SCI31_I: return R[2];
			case RNFRSCal::SCI41_I: return R[3];
			default: mnd::unreachable();
		}
	}
	std::unique_ptr<phy::Nucleus> SecondaryIon() const;
	void GetFrom(const mnd::RunsheetState& );

    BeamInfo() = default;
	virtual ~BeamInfo() = default;
	ClassDef(BeamInfo, 1);
};
ADD_JSON_TYPE_RESOLUTION(BeamInfo, 3)

struct RNFRSHit {
	static constexpr i32 S2_BT_TRACKING_INCLUDE_SCI21_MASK            = 0x01;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC21_MASK            = 0x02;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC22_MASK            = 0x04;
	static constexpr i32 S2_BT_TRACKING_INCLUDE_TPC23_MASK            = 0x08;

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
        f64 AoQ  = NAN;
        f64 Q    = NAN;
		f64 x0   = NAN; // at nominal z = 0 ; FRS standard coordinates
		f64 y0   = NAN; // at nominal z = 0 ; FRS standard coordinates
		f64 ax   = NAN;
		f64 ay   = NAN;
		f64 beta = NAN;

		i64 code = 0; // some metadata fed from the processor.
		
		Id() = default;

		inline void Clean() noexcept {
			AoQ=Q=x0=y0=ax=ay=beta = NAN;
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
	}

	virtual ~RNFRSHit() = default;
	ClassDef(RNFRSHit, 1);
};

struct TFRSHitCont : TContainer<RNFRSHit> {
	TH2D *h2_track_x;
	TH2D *h2_track_y;
	TH2D *h2_target_xy;

	TH2D *h2_s2_q;
	TH2D *h2_s2_qc;

	TH2D *h2_s3_q;
	TH2D *h2_s3_qc;

	TH1D *h1_s3_beta;
	TH2D *h2_s3_id;

	/* Next fields just get forwarded by the processor, from the cal step. */
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param{}; // gets FW'ed from cal step.
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_param{}; // gets FW'ed from cal step.
    TrigParam *trig_param;                                    // gets FW'ed from cal step.
	std::string *setupName;                                   // gets FW'ed from cal step.

    FRSToFParam *sTof;
	BeamInfo *binfo;
	
	TFRSHitCont();

	void Setup() override;
};

namespace mnd::geom { struct Line3D; }
mnd::geom::Line3D RNTrackToLine3D(const RNFRSHit::Id& );
