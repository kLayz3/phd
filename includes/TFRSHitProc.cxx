
#include "TFRSCalCont.h"
#include "TFRSHitCont.h"
#include "TFRSHitProc.h"
#include "util/PolyFitter.h"
#include "util/json_struct_def.hh"
#include <cmath>

#include "util/Geometry.h"

constexpr double TWO_PI = 6.283185307179586;

thread_local mnd::geom::Point2D g_upstream_vertex {};

#ifdef _MSC_VER
#	include <intrin.h>
#endif

/* Count bits set in an integer. */
[[maybe_unused]]
static int count(int v) noexcept {
#ifndef _MSC_VER
	return __builtin_popcount(v);
#else
	return __popcnt(v);
#endif
}

/* Return a number in interval: [0,1> */
static double uniform() noexcept {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1ULL);
}

/* Draw a random point (x,y) inside a circle of radius `R` */
static std::array<double, 2> uniform_circle(double R) noexcept {
	double t = uniform();
	double r = R * sqrt(t);
	double phi = TWO_PI * uniform();
	return { r*cos(phi), r*sin(phi) };
}

TFRSHitProc::TFRSHitProc (
	TFRSHitCont& out,
	const TFRSCalCont& in, 
	int s2_bt_mask
) : TFRSHitProc::Base(out, in),
	s2_bt_tracking_mask(s2_bt_mask)
{
	x.reserve(8); y.reserve(8);
	zx.reserve(8); zy.reserve(8);
	
	extern double g_expert_target_z; // extern'ed from: `includes/TFOOTHitCont.cxx`
	z0        = g_expert_target_z;
	tar_width = out.sTar->width;

	out.h2_target_xy->SetTitle(Form("%s@z0=%.1f", out.h2_target_xy->GetTitle(), z0));
}

void TFRSHitProc::ProcessEntry() noexcept {
	RNFRSHit& out = (this->out).inner();
	out.Clean();

	out.cal = std::get<0>(this->in).inner(); // RNFRSCal& operator=(RNFRSCal& )
	ProcessS2BT();
}

void TFRSHitProc::ProcessS2BT() noexcept {
	x.clear(); y.clear();;
	zx.clear(); zy.clear();;

	const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	
	/* Local references. */
	RNFRSHit::Id& bt = out.inner().s2_bt;
	double& xT = out.inner().xT; 
	double& yT = out.inner().yT;

	/* This is a mask to say that for incoming track, a (0,0) point directly on the target 
	 * also gets included. Just for debugging / sanity checks. Do not use for real data. */
	if(s2_bt_tracking_mask & RNFRSHit::S2_BT_TRACKING_INCLUDE_POINTLIKE_TARGET_MASK) {
		auto [x0,y0] = uniform_circle(tar_width);
		x.push_back(x0);
		y.push_back(y0);
		zx.push_back(z0);
		zy.push_back(z0);
		bt.code |= RNFRSHit::S2_BT_TRACKING_INCLUDE_POINTLIKE_TARGET_MASK;
	}

	if(s2_bt_tracking_mask & RNFRSHit::S2_BT_TRACKING_INCLUDE_SCI21_MASK) {
		static const double sci21_z = cal.sci_param->at(0).z0;
		const std::vector<RNSciCal::Measurement>& hits = in.sci[0].hits; 
		if(hits.size() == 1) { // otherwise can't resolve.
			x.push_back( hits[0].x );
			zx.push_back( sci21_z );
			bt.code |= RNFRSHit::S2_BT_TRACKING_INCLUDE_SCI21_MASK;
		}
	}

#define TRY_INCLUDE_TPC_INTO_BT_TRACKING(LABEL, INDEX) \
	if(s2_bt_tracking_mask & RNFRSHit::S2_BT_TRACKING_INCLUDE_TPC##LABEL##_MASK) { \
		const double ztpc = cal.tpc_param->at(INDEX).z0; \
		const RNTPCCal& tpc = in.tpc[INDEX]; \
		const double xtpc = tpc.X0(); \
		const double ytpc = tpc.Y0(); \
		if( std::isfinite(xtpc) and std::isfinite(ytpc) ) { \
			x.push_back( xtpc ); \
			y.push_back( ytpc ); \
			zx.push_back( ztpc ), zy.push_back( ztpc ); \
			\
			bt.code |= RNFRSHit::S2_BT_TRACKING_INCLUDE_TPC##LABEL##_MASK; \
		} \
	} \
	EMPTY_MACRO__(INDEX)
	
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(21, 0)
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(22, 1)
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(23, 2)

	if(x.size() >= 2) {
		PolyFit<1>(zx, x, this->fit_result);
		bt.x0 = fit_result[0]; 
		bt.ax = fit_result[1];
		xT = poly::Eval(z0, fit_result);
	}
	if(y.size() >= 2) {
		PolyFit<1>(zx, y, this->fit_result);
		bt.y0 = fit_result[0];
		bt.ay = fit_result[1]; 
		yT = poly::Eval(z0, fit_result);
	}
	
	out.h2_target_xy->Fill(xT, yT);
	
	g_upstream_vertex = { .x = xT, .y = yT };
}
