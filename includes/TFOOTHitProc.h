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
	static constexpr double DEFAULT_COST_R = 100.0; // default cost per mm^2 difference
	static constexpr double DEFAULT_COST_Q = 10.0;  // default cost per charge^2 difference (not sure!)
	static constexpr double DEFAULT_COST_T =  4.0;  // default cost per mm^2 of upstream @target difference
	static constexpr double DEFAULT_COST_P = 1e10;  // default cost if next layer missed
	/* For 12C the average sqrt(kq) cost is ~0.3 or so, so variance is ~0.1 or so.
	 * sqrt(kr) is anything between 1-5mm (worst case, probably I fucked up alignment. 
	 * Should be ~1mm on a good day. */

	// For k==3 degrees of freedom: Chi^2;
	constexpr static double DEFAULT_MAX_COST = 15;

	/* Maximum cost to allow a point to branch the Kalman cannot be constant,
	 * as it depends on the layer number:
	 * Layer [1] => only `kq` cost
	 * Layer [2] => `kq` + `kt` cost
	 * Layer [3] => `kq` + `kt` + `kr` cost
	 * Layer [4] => `kq` + `kt` + `kr` (3-point-fit) cost */
	struct MaxCost {
		MaxCost() = default;
		MaxCost(double x) : overall_cost(x) {}
		inline double at(size_t n) const noexcept {
			switch(n) {
				case 0: return overall_cost * 4.0;
				case 1: return overall_cost * 3.0;
				case 2: return overall_cost * 1.3;
				case 3: return overall_cost * 1.0;
				default:
					__builtin_unreachable();

			};
		}
		inline double operator*() const noexcept { return overall_cost; }

	private:
		double overall_cost = TrackCost::DEFAULT_MAX_COST;
	};
	enum F { KR, KQ, KP, KT };

	/* Exposing setters and getters because the sideffect is bookkeeping
	 * the sum on the fly, which raw reference accesses would invalidate :-) */
	inline double kr() const noexcept { return kr_; }
	inline double kq() const noexcept { return kq_; }
	inline double kp() const noexcept { return kp_; }
	inline double kt() const noexcept { return kt_; }

	/* Set the value of individual cost component, and update the total sum.
	 * Handles NAN's and is idempotent. */
	template<enum F o>
	void set(double v) noexcept {
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
		
		if(std::isfinite(v))
			*sum_ += v;
	}

	inline double sum() const noexcept { return sum_ ? *sum_ : std::numeric_limits<double>::infinity(); }
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
	using FTrackOnline = Track<N_PAIRS, RNFOOTPair>;

	constexpr static double CLUSTER_SIZE_ONE_Q_CUTOFF = 1.5; // when cluster size == 1 doesn't make sense anymore.
	constexpr static double TARGET_Z = 0.0; // by convention. In Kalman coordinates, place target nominally at 0.0
											// Will be shifted back to "real" FRS coordinates later.
	constexpr static u32 MAX_CANDIDATES = 10;
	using DAG = DirectedAGraph<u16, N_PAIRS>;

	TFOOTHitProc(TFOOTHitCont& , BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (const,&) ), 
		double = TrackCost::DEFAULT_MAX_COST,
		const std::array<double,4>& = {NAN, NAN, NAN, NAN}, /* cost coefficients: {Cr, Cq, Ct, Cp} */
		bool  = false, /* require_valid_upstream_track */
		Verbosity = Verbosity::SILENT);
	TFOOTHitProc() = default;

	double q_tolerance;
	TrackCost::MaxCost max_cost;
	static Verbosity v;	

	double Cr = TrackCost::DEFAULT_COST_R; // cost per mm^2 difference
	double Cq = TrackCost::DEFAULT_COST_Q; // cost per charge^2 difference
	double Ct = TrackCost::DEFAULT_COST_T; // cost per mm^2 of upstream @target difference
	double Cp = TrackCost::DEFAULT_COST_P; // cost if next layer missed

	double kr(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 
	double kq(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 
	std::pair<double, double> kt_kp(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 

	void ProcessEntry() noexcept;

private:
	void ProcessPair(const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>&, i32) noexcept;
	void ConstructObviousTracks() noexcept;
	void ConstructDAG() noexcept;
	void AnalyseDAG() noexcept;
	void PreProcess() noexcept;
	void PostProcess() noexcept;

	bool requires_valid_upstream_track;

	std::array<double, N_PAIRS> pair_z;
	mnd::geom::Rectangle2D target_xy;
	mnd::geom::Point2D upstream_hit_loc;
	std::vector<mnd::geom::Line3D> lines{};

	DAG dag;
	
	std::array<FHitMatrix, N_PAIRS> hm; // hit matrices
	std::array<Eigen::Matrix2d, N_PAIRS> A_inv; // A_inv[i] == hm[i]::A inverted
	std::array<Eigen::Vector2d, N_PAIRS> refl; // +-1 based on the on each of the `orientation` params  
	void SetConversionMatrices(int , const FOOTParam& , const FOOTParam& );

	std::array <
		std::pair<double, DAG::Index>, 50
	> path_specific_candidates_buf;

	FTrackOnline GetPrelimTrackFromPath(const DAG::DAGPath& ) const noexcept;
	void PoisonEntriesFromHMs(const DAG::DAGPath&) noexcept; 
};

