#pragma once

#include "Eigen/Core"
#include "monad/monad.hxx"

#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <cmath>
#include <numeric>

#include "util/Geometry.h"


#define GEN_ARG_TYPE_FOOT(z, n, data) \
	BOOST_PP_TUPLE_ELEM(2, 0, data) TFOOTCalCont BOOST_PP_TUPLE_ELEM(2, 1, data)

struct TFOOTHitProc : TProcessor <
	TFOOTHitCont
	( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )
> {
	using Base = TProcessor<TFOOTHitCont( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )>;
	using FHitMatrix = HitMatrix<RNFOOTPair>;
	using FTrack = Track<RNFOOTPair>;

	constexpr static double DETECTOR_MIDPOINT = 319.5;
	constexpr static double STRIP_TO_MM       = 0.150;

	constexpr static double DEF_TOLERANCE = 0.3;
	constexpr static u32 N_PAIRS = TFOOTHitCont::N_PAIRS;

	TFOOTHitProc(TFOOTHitCont& , BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (const,&) ), double t = DEF_TOLERANCE );
	TFOOTHitProc() = default;

	void ProcessEntry() noexcept;

	// If the dependence E(Z) = A * Z^a, then:
	// f = 1/A, c = 1/a <=> Z(E) = (f * E)^a 
	struct e_to_z_t {
		double f, c;
		inline double operator()(double e) const noexcept {
			return std::pow(f*e, c);
		}
		void Init(const FOOTParam& );
	}; 
	std::array<e_to_z_t, N_FOOT_DETECTORS> e_to_z;

	double Cr, CQ, fT, cP;
	double kr(const FTrack& , const FHitMatrix::Entry& , int ); 
	double kQ(const FTrack& , const FHitMatrix::Entry& , int ); 
	double kt_kp(const FTrack& , const FHitMatrix::Entry& , int ); 
	double k(const FTrack& , const FHitMatrix::Entry& , int );

private:
	void ProcessPair(const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>&, i32) noexcept;
	void ProcessTracks() noexcept;

	std::array<double, N_PAIRS+1> pair_z; // last layer is a no-op
	
	double target_z;
	mnd::geom::Rectangle2D target_area;

	
	std::array<FHitMatrix, N_PAIRS> hm;
};

