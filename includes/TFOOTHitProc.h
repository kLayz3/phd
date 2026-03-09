#pragma once

#include "monad/monad.hxx"

#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#define GEN_ARG_TYPE_FOOT(z, n, data) \
	BOOST_PP_TUPLE_ELEM(2, 0, data) TFOOTCalCont BOOST_PP_TUPLE_ELEM(2, 1, data)

struct TFOOTHitProc : TProcessor <
	TFOOTHitCont
	( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )
> {
	using Base = TProcessor<TFOOTHitCont( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )>;
	constexpr static double DEF_TOLERANCE = 0.3;

	TFOOTHitProc(TFOOTHitCont& , BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (const,&) ), double t = DEF_TOLERANCE );
	TFOOTHitProc() = default;

	void ProcessEntry() noexcept;

private:
	void ProcessSingle(const TFOOTCalCont& ) noexcept;
	bool IsCompatible(const double, const double) noexcept; 

	double e_diff_tolerance;

	struct HitCandidate {
		double e, pos;
		HitCandidate(double e_, double pos_) : e(e_), pos(pos_) {};
		HitCandidate() = default;
	};
	
	struct HitsBuffer {
		constexpr static double DIFF_ALLOWED = 0.1;

		int nx = -1; /* CFOOT index for the x-measuring FOOT */
		int ny = -1; /* CFOOT index for the y-measuring FOOT */
		double zx = NAN; /* z- value passed in from the parameter file. */
		double zy = NAN; /* z- value passed in from the parameter file. */
		std::vector<HitCandidate> xs, ys;
		inline void Clear() noexcept { nx = -1; ny = -1; zx = NAN; zy = NAN; xs.clear(), ys.clear(); }
		inline bool IsFilled() const noexcept {
			return nx > -1 and ny > -1 and
				nx != ny and std::isfinite(zx) and std::isfinite(zy) and
				xs.size() > 0 and ys.size() > 0;
		}
		/* Just checks that both xz and zy correspond to same value. */
		inline bool IsValid() const noexcept {
			if(!IsFilled()) return false;

			if( std::abs(zx - zy) > DIFF_ALLOWED ) 
				ERROR("Validly filled FOOT pair, but z-distance mismatched? "
					"nx = %d, zx = %.3f; ny = %d, zy = %.3f\n", nx,zx,ny,zy);
			return true;
		}
		inline double Z() const noexcept {
			return (zx + zy) / 2;
		}
	};
	
	void ConstructHits(HitsBuffer&, int) noexcept;

	HitsBuffer buf[ TFOOTHitCont::N_PAIRS ];
};

