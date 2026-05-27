#pragma once

#include "monad/monad.hxx"
#include "Eigen/Core"

#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <cmath>
#include <limits>
#include <numeric>

#include "util/Geometry.h"
#include "util/DirectedAGraph.hxx"
#include "util/Verbosity.hxx"

#define GEN_ARG_TYPE_FOOT(z, n, data) \
	BOOST_PP_TUPLE_ELEM(2, 0, data) TFOOTCalCont BOOST_PP_TUPLE_ELEM(2, 1, data)

/* Container to keep track of the cost that the track candidate acquires over the
 * CKF algorithm. Real value only comes from the sum, but keeping individual components
 * checked is used to e.g. normalize different coefficients... */
struct TrackCost {
	static constexpr double DEFAULT_COST_R = 1.0;  // default cost per mm^2 difference
	static constexpr double DEFAULT_COST_Q = 5.0;  // default cost per charge^2 difference
	static constexpr double DEFAULT_COST_T = 1e10; // default cost if target missed
	static constexpr double DEFAULT_COST_P = 1e2;  // default cost if next layer missed

	enum F { KR, KQ, KP, KT };

	/* Exposing setters and getters because the sideffect is bookkeeping
	 * the sum on the fly, which raw reference accesses would invalidate :-) */
	inline double kr() const noexcept { return kr_; }
	inline double kq() const noexcept { return kq_; }
	inline double kp() const noexcept { return kp_; }
	inline double kt() const noexcept { return kt_; }
	
	template<enum F o>
	void set(double v) noexcept {
#ifdef MND_HITMATRIX_DO_BOUNDS_CHECK
		assert(std::isfinite(v) && "Must pass a finite value here.");
#endif
		if(!sum_) sum_ = 0;

		if constexpr(o == KR) {
			if(std::isfinite(kr_)) *sum_ -= kr_;
			kr_ = v;
		}
		if constexpr(o == KQ) {
			if(std::isfinite(kq_)) *sum_ -= kq_;
			kq_ = v;
		}
		if constexpr(o == KP) {
			if(std::isfinite(kp_)) *sum_ -= kp_;
			kp_ = v;
		}
		if constexpr(o == KT) {
			if(std::isfinite(kt_)) *sum_ -= kt_;
			kt_ = v;
		}
		*sum_ += v;
	}

	double sum() const noexcept { return sum_ ? *sum_ : std::numeric_limits<double>::infinity(); }
	friend std::ostream& operator<<(std::ostream&, const TrackCost&);

private:
	double kr_ = NAN, kq_ = NAN, kp_ = NAN, kt_ = NAN;
	std::optional<double> sum_ {};
};

struct TFOOTHitProc : TProcessor <
	TFOOTHitCont
	( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )
> {
	using Base = TProcessor<TFOOTHitCont( BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (,) ) )>;
	constexpr static u32 N_PAIRS = TFOOTHitCont::N_PAIRS;

	using FHitMatrix = HitMatrix<RNFOOTPair>;
	using FTrack = Track<N_PAIRS, RNFOOTPair>;

	constexpr static int N_STRIPS = TFOOTCalCont::N_STRIPS;
	constexpr static double DETECTOR_MIDPOINT = static_cast<double>(N_STRIPS - 1) / 2; // strip units
	constexpr static double DETECTOR_SIZE     = 96.0; // [ mm ]
	constexpr static double STRIP_TO_MM       = DETECTOR_SIZE / N_STRIPS; // 0.150 mm/strip
	constexpr static double MM_TO_STRIP       = 1.0 / STRIP_TO_MM; // 6.67 strip/mm

	constexpr static double CLUSTER_SIZE_ONE_Q_CUTOFF = 1.5; // when cluster size == 1 doesn't make sense anymore.
	constexpr static double DEFAULT_MAX_Q_TOLERANCE = 0.8;
	constexpr static double DEFAULT_MAX_COST = 100;
	constexpr static double TARGET_Z = 0.0; // convention; safe to hard-bake

	using DAG = DirectedAGraph<u16, N_PAIRS>;

	TFOOTHitProc(TFOOTHitCont& , BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (const,&) ), 
		double = DEFAULT_MAX_Q_TOLERANCE,
		double = DEFAULT_MAX_COST,
		Verbosity = Verbosity::SILENT);
	TFOOTHitProc() = default;

	// If the dependence E(Q) = A * Q^a, then:
	// f = 1/A, c = 1/a <=> Q(E) = (f * E)^a 
	struct e_to_q_t {
		double f, c;
		inline double operator()(double e) const noexcept {
			return std::pow(f*e, c);
		}
		void Init(const FOOTParam& );
	}; 
	std::array<e_to_q_t, N_FOOT_DETECTORS> e_to_q;
	
	double q_tolerance;
	double max_cost;
	static Verbosity v;	

	double Cr = TrackCost::DEFAULT_COST_R;
	double Cq = TrackCost::DEFAULT_COST_Q; 
	double Ct = TrackCost::DEFAULT_COST_T;
	double Cp = TrackCost::DEFAULT_COST_P;

	double kr(const FTrack& , const FHitMatrix::Entry& , u32 ) const; 
	double kq(const FTrack& , const FHitMatrix::Entry& , u32 ) const; 
	std::pair<double, double> kt_kp(const FTrack& , const FHitMatrix::Entry& , u32 ) const; 

	void ProcessEntry() noexcept;

private:
	void ProcessPair(const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>&, i32) noexcept;
	void ConstructObviousTracks() noexcept;
	void ConstructDAG() noexcept;

	std::array<double, N_PAIRS> pair_z;
	mnd::geom::Rectangle2D target_xy;

	DAG dag;
	
	std::array<FHitMatrix, N_PAIRS> hm; // hit matrices
	std::array<Eigen::Matrix2d, N_PAIRS> A_inv; // A_inv[i] == hm[i]::A inverted
	std::array<Eigen::Vector2d, N_PAIRS> refl; // +-1 based on the on each of the `orientation` params  
	void SetConversionMatrices(int , const FOOTParam& , const FOOTParam& );

	FTrack GetPrelimTrackFromPath(const DAG::DAGPath& ) const;
};

