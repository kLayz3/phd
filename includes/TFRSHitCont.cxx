#include "TFRSHitCont.h"
#include "TFRSCalCont.h"

#include "util/Geometry.h"
#include "util/MPhysics.h"
#include "util/PolyFitter.h"
#include "ExperimentAssumptions.hh"
#include "util/RunsheetParser.h"

using nlohmann::json;

static std::array<TPCParam, RNFRSCal::N_VALID_TPC> _tpc_param {};
static std::array<SCIParam, RNFRSCal::N_VALID_SCI> _sci_param {};

double FRSToFSingle::Beta(double dt) const noexcept {
	return 1.0 / poly::Eval(dt, par);
}

double FRSToFSingle::Beta(const RNFRSCal& cal) const noexcept {
	const RNSciCal& start = cal.sci[ combo[0] ];
	const RNSciCal& end = cal.sci[ combo[1] ];
	/* ^^^ both acceses shall not be UB.
	 * Validity checked during `TFRSHitProc::TFRSHitProc(..)` ctor */

	if(start.hits.empty() || end.hits.empty())
		return NAN;
	
	const RNSciCal::Measurement& m0 = start.hits.front();
	const RNSciCal::Measurement& m1 = end.hits.front();
	
	return this->Beta(m1.t - m0.t);
}
FRSToFSingle const* FRSToFParam::Get(u32 start, u32 stop) const noexcept {
	for(const auto& attempt : this->ToF) {
		const auto& combo = attempt.combo;
		if(start == combo[0] && stop == combo[1]) {
			return &attempt;
		}
	}
	return nullptr;
}

bool FRSToFSingle::IsValid() const noexcept {
	return combo[0] < combo[1]
		&& combo[1] < (u32)RNFRSCal::N_VALID_SCI
		&& par[1] > 0;
}

void BeamInfo::GetFrom(const mnd::RunsheetState& runsheet_row) {
	A0 = runsheet_row.secondary.A;
	Z0 = runsheet_row.secondary.Z;
	Q0 = runsheet_row.secondary.Q();

	R[0] = runsheet_row.brho[mnd::TA_S1];
	R[1] = runsheet_row.brho[mnd::S1_S2];
	R[2] = runsheet_row.brho[mnd::S2_S3];
	R[3] = runsheet_row.brho[mnd::S3_S4];
}

std::unique_ptr<phy::Nucleus> BeamInfo::SecondaryIon() const {
	mnd::Option<u16> ne = mnd::None;
	if(Z0 >= Q0)
		ne = mnd::Some{(u16)(Z0 - Q0)};
	return std::unique_ptr<phy::Nucleus>(
		new phy::Nucleus {
		  .A = this->A0,
		  .Z = this->Z0,
		  .N_electrons = ne
		}
	);
}

std::string RNFRSHit::DecodeS2() const noexcept {
	std::stringstream info;
	info << "S2 BT tracking: ";
	if(s2_bt.code & S2_BT_TRACKING_INCLUDE_SCI21_MASK) {
		info << "SCI21 ";
	}
	if(s2_bt.code & S2_BT_TRACKING_INCLUDE_TPC21_MASK) {
		info << "TPC21 ";
	}
	if(s2_bt.code & S2_BT_TRACKING_INCLUDE_TPC22_MASK) {
		info << "TPC22 ";
	}
	if(s2_bt.code & S2_BT_TRACKING_INCLUDE_TPC23_MASK) {
		info << "TPC23 ";
	}
	return info.str();
}
std::string RNFRSHit::DecodeS3() const noexcept {
	return "Work in progress! 🚧";
}

TFRSHitCont::TFRSHitCont() : TContainer("FRS") {}

void TFRSHitCont::Setup() {
	h2_track_x   = RegisterObject<TH2D>("s2_upstream_track_x", "S2 TPC Tracking;z[mm];x[mm]", 400, 0, 4200, 200, -100, 100);
	h2_track_y   = RegisterObject<TH2D>("s2_upstream_track_y", "S2 TPC Tracking;z[mm];y[mm]", 400, 0, 4200, 200, -100, 100);
	h2_target_xy = RegisterObject<TH2D>("s2_target_xy", "S2 Upstream Target Hit;x[mm];y[mm]", 200, -20, 20, 200, -20, 20);

	h2_s2_q      = RegisterObject<TH2D>("s2_q", "S2 Particle Charge;Q[before S2 target];Q[after S2 target]",
		200, 0, 10, 200, 0, 10);
	h2_s2_qc     = RegisterObject<TH2D>("s2_qc", "S2 Particle Charge after #beta-corr (both measurements);Q[before S2 target];Q[after S2 target]",
		200, 0, 10, 200, 0, 10);
	
	h2_s3_q      = RegisterObject<TH2D>("s3_q", "S2 Particle Charge;Q[before S2 target];Q[S3]",
		200, 0, 10, 200, 0, 10);
	h2_s3_qc     = RegisterObject<TH2D>("s3_qc", "S2 Particle Charge after #beta-corr (both measurements);Q[before S2 target];Q[S3]",
		200, 0, 10, 200, 0, 10);

	h1_s3_beta   = RegisterObject<TH1D>("s3_b", "S3 beta velocity;#beta = v/c", 400, 0, 1);
	h2_s3_id     = RegisterObject<TH2D>("s3_id", "S3 Particle ID;AoQ;Q[charge]",
		200, 1.0, mnd::assume::primary.AoQ() + 0.2,
		200, 0.0, mnd::assume::primary.Z + 2
	);
	
	tpc_param  = RegisterObject<std::array<TPCParam, RNFRSCal::N_VALID_TPC>>("tpc_parameters", {});
	sci_param  = RegisterObject<std::array<SCIParam, RNFRSCal::N_VALID_SCI>>("sci_parameters", {});
	trig_param = RegisterObject<TrigParam>("trigger_map", TrigParam{});
	setupName  = RegisterObject<std::string>("setup_file", {});

	sTof  = RegisterObject<FRSToFParam>("tof_par", mnd::noop_fn<FRSToFParam>(), FRSToFParam{}); // Will be filled out by `TFRSHitProc::TFRSHitProc(..) ctor`
	binfo = RegisterObject<BeamInfo>("beam_info" , mnd::noop_fn<BeamInfo>(),    BeamInfo{});    // Will be filled out by `TFRSHitProc::TFRSHitProc(..) ctor`
}

mnd::geom::Line3D RNTrackToLine3D(const RNFRSHit::Id& t) { return { t.x0, t.ax, t.y0, t.ay }; }

ClassImp(FRSToFSingle);
ClassImp(FRSToFParam);
ClassImp(BeamInfo);
ClassImp(RNFRSHit);
ClassImp(RNFRSHit::Id);
