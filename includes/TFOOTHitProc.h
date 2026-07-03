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
	static constexpr double DEFAULT_COST_Q =  12.0; // default cost per charge^2 difference (TODO: not sure)
	static constexpr double DEFAULT_COST_T =   6.0; // default cost per mm^2 of upstream @target difference
	static constexpr double DEFAULT_COST_P = 100.0; // default cost if next layer missed
	/* For 12C the average sqrt(kq) cost is ~0.3 or so, so variance is ~0.1 or so.
	 * sqrt(kr) is anything between 1-5mm (worst case, probably I messed up alignment. 
	 * Should be ~1mm on a good day. */

	// For k==3 degrees of freedom: total Chi^2 of about 12 is 99% confidence. 
	// In ideal world... but both `kt` and `kq` can dance like crazy. */
	constexpr static double DEFAULT_MAX_CANDIDATE_COST = 1'000;
	constexpr static double DEFAULT_MAX_FINAL_COST = 22;

	constexpr static double NIL_VALUE = NAN;
	
	enum F { KR, KQ, KP, KT };

	static double Cr;
	static double Cq;
	static double Ct;
	static double Cp;
	static double max_cost;
	static double max_cost_final_track;

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
	template<enum F o>
	bool is_set() const noexcept { 
		if constexpr(o == KR) {
			return std::isfinite(kr_);
		}
		else if constexpr(o == KQ) {
			return std::isfinite(kq_);
		}
		else if constexpr(o == KP) {
			return std::isfinite(kp_);
		}
		else if constexpr(o == KT) {
			return std::isfinite(kt_);
		}
	}

	inline double sum() const noexcept { return sum_ ? *sum_ : std::numeric_limits<double>::infinity(); }
	inline void reset() noexcept { sum_.reset(); }

	friend std::ostream& operator<<(std::ostream&, const TrackCost&);
	friend bool operator>(const TrackCost& , double );

	TrackCost() = delete;
	TrackCost(u32 n) : test_track_size(n) {};
	
	u32 test_track_size;

private:
	double kr_ = NAN, kq_ = NAN,  kt_ = NAN, kp_ = NAN;
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

	constexpr static float CLUSTER_SIZE_ONE_Q_CUTOFF = 1.7; // when cluster size == 1 doesn't make sense anymore.
	constexpr static double TARGET_Z = 0.0; // by convention. In Kalman coordinates, place target nominally at 0.0
											// Will be shifted back to "real" FRS coordinates later.
	constexpr static u32 MAX_CANDIDATES = 10; // To how many paths can a node branch to (at most)
	using DAG = DirectedAGraph<u16, N_PAIRS>;

	TFOOTHitProc(TFOOTHitCont& , BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_TYPE_FOOT, (const,&) ), 
		double = TrackCost::DEFAULT_MAX_CANDIDATE_COST,
		double = TrackCost::DEFAULT_MAX_FINAL_COST,
		const std::array<double,4>& = {NAN, NAN, NAN, NAN}, // cost coefficients: {Cr, Cq, Ct, Cp}
		bool  = false,                                      // requires_valid_upstream_track
		Verbosity = Verbosity::SILENT);
	TFOOTHitProc() = default;

	static Verbosity v;	
	static bool requires_valid_upstream_track;

	double kr(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 
	double kq(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 
	std::pair<double, double> kt_kp(const FTrackOnline& , const FHitMatrix::Entry& , u32 ) const; 

	void ProcessEntry() noexcept;

	using CandidatesBuffer = std::array <
		std::pair<double, DAG::Index>, 50
	>;

private:
	void ProcessPair(const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>&, i32) noexcept;
	void ConstructObviousTracks() noexcept;
	void ConstructDAG() noexcept;
	void AnalyseDAG() noexcept;
	void PreProcess() noexcept;
	void PostProcess() noexcept;

	std::array<double, N_PAIRS> pair_z;
	mnd::geom::Rectangle2D target_xy;
	mnd::geom::Point2D upstream_hit_loc;
	std::vector<mnd::geom::Line3D> lines{};

	DAG dag;
	
	std::array<FHitMatrix, N_PAIRS> hm; // hit matrices
	std::array<Eigen::Matrix2d, N_PAIRS> A_inv; // A_inv[i] == hm[i]::A inverted
	std::array<Eigen::Vector2d, N_PAIRS> refl; // +-1 based on the on each of the `orientation` params  
	void SetConversionMatrices(int , const FOOTParam& , const FOOTParam& );

	CandidatesBuffer path_specific_candidates_buf;

	FTrackOnline GetPrelimTrackFromPath(const DAG::DAGPath& ) const noexcept;
	void PoisonEntriesFromHMs(const DAG::DAGPath&) noexcept; 
};

