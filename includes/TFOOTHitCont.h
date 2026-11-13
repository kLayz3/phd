#pragma once

#include "monad/monad.hxx"
#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "TFOOTCalCont.h"

class TH2I;

/* Heavy debugging features compiled in. Can be toggled together automatically.
 * This is either for sanity debugging or for:
 *  - quick/hit/foot_track_analysis.cc ,
 *  - quick/hit/foot_regain_match.cc
 */
//#define MND_FOOTTRACK_DEBUG

/* In case the heavy debug build is enabled, also feature it in here. */
#if defined(MND_DEBUG_ENABLED) && !defined(MND_FOOTTRACK_DEBUG)
#	define MND_FOOTTRACK_DEBUG
#endif

#ifdef MND_FOOTTRACK_DEBUG
#	ifndef MND_HITMATRIX_DO_BOUNDS_CHECK
#		define MND_HITMATRIX_DO_BOUNDS_CHECK
#	endif
#   pragma message(MND_FILE_LINE_STR " => FOOT Hit Step Debugging Features Enabled!")
#endif

#include "util/HitMatrix.hxx"
#include "util/FTrack.h"

namespace mnd::geom { struct Line3D; }

/* Represent the 'charge' measurement of each layer. */
struct FOOTQ {
	using ClusterType = RNFOOTCluster::ClusterType;

	float q; /* Calibrated charge value [e] */
	u32 fCM = 0; /* Cluster multiplicity. */

	/* Few fields taken from RNFOOTCluster.. */
#ifdef MND_FOOTTRACK_DEBUG
	ClusterType fCT; /* Cluster type. */
	float delta; /* Delta coefficient of the cluster. */
	float e0; /* Uncorrected cluster energy. */
#endif
	
	FOOTQ() = default;
	FOOTQ(float q_, u32 fCM_
	#ifdef MND_FOOTTRACK_DEBUG
		, ClusterType fCT_, float d_, float e0_
	#endif
	) :
		q(q_), fCM(fCM_)
	#ifdef MND_FOOTTRACK_DEBUG
		, fCT(fCT_), delta(d_), e0(e0_)
	#endif
	{}

	virtual ~FOOTQ() = default;
	ClassDef(FOOTQ, 1);
};
#ifndef MND_FOOTTRACK_DEBUG
static_assert(sizeof(FOOTQ) == 16, "Alignment obsession, friends.");
#endif

struct FOOTHit {
	using ClusterType = FOOTQ::ClusterType;

	FOOTQ Q;
	double m; // measurement [ mm ], in the writeup corresponds to either `u` or `v` .
	          // it is *only* converted to mm. There is no FOOT offset added here.
	
	FOOTHit() = default;
	FOOTHit(double q_, u32 fCM_
	#ifdef MND_FOOTTRACK_DEBUG
		, ClusterType fCT_, float d_, float e0_
	#endif
		, double m_) :
		Q(q_, fCM_
	#ifdef MND_FOOTTRACK_DEBUG
			, fCT_, d_, e0_
	#endif
		), m(m_) {}

	virtual ~FOOTHit() = default;
	ClassDef(FOOTHit, 1);
};
#ifndef MND_FOOTTRACK_DEBUG
static_assert(sizeof(FOOTHit) == 32, "Alignment obsession, friends.");
#endif
std::ostream& operator<<(std::ostream& , const FOOTHit& ) noexcept;

struct RNFOOTPair {
	std::vector<FOOTHit> x;
	std::vector<FOOTHit> y;

    /* Use the lowest bit of the mantissa to indicate if the
     * `x` or the `y` foot are missing. This changes the precision by ~1e-14 or so
     * such that it doesn't really matter. Sane person would implement a separate flag, but I'm an
     * OCD masochist who likes alignment very much. Not safe for kids. :) */
	f32 zx, zy;

	RNFOOTPair() = default;

    constexpr static u32 _DATA_PRESENT_BITINDEX = 0;
    constexpr static u64 _DATA_PRESENT_BITMASK  = (1ULL << _DATA_PRESENT_BITINDEX);

    bool HasDataX() const noexcept;
    bool HasDataY() const noexcept;
    bool HasData() const noexcept;

	inline constexpr double GetZx() const noexcept { return static_cast<double>(zx); };
	inline constexpr double GetZy() const noexcept { return static_cast<double>(zy); };

	/* No need to clean z's. They are static for a sequence of files. */
	inline void Clean() noexcept { x.clear(); y.clear(); }
	virtual ~RNFOOTPair() = default;
	ClassDef(RNFOOTPair, 1);
};
static_assert(sizeof(RNFOOTPair) == 64, "Alignment obsession, friends.");

struct RNFOOTTrack {
	/* Lay them out in memory just like the Normal repr of Line3D expects. */
	double x0, ax;
	double y0, ay;
	double Q;

	double score;
	std::size_t n = 0; // number of collected pairs, size_t anyway chosen as it will be aligned to 8-byte

	inline bool operator<(const RNFOOTTrack& rhs) const noexcept {
		return Q > rhs.Q; // descending charge.
	}

	/* Comparison and deref operators required for Vertexing API in Geometry.h */
	friend bool operator==(const RNFOOTTrack& , const RNFOOTTrack &) noexcept;
	inline constexpr mnd::geom::Line3D operator*() const noexcept { return {x0,ax,y0,ay}; };
	inline constexpr mnd::geom::Line3D::Normal Normal() const noexcept { return {x0,ax,y0,ay}; }

#ifdef MND_FOOTTRACK_DEBUG
	static constexpr u32 N_PAIRS = N_FOOT_DETECTORS / 2;
	
	std::array<double, N_PAIRS> _u;
	std::array<double, N_PAIRS> _v;
	std::array<double, N_PAIRS> _zx;
	std::array<double, N_PAIRS> _zy;
	std::array<double, N_PAIRS> _qx;
	std::array<double, N_PAIRS> _qy;

	std::array<double, N_PAIRS> _e0_x; // accumulated cluster energy
	std::array<double, N_PAIRS> _delta_x;
	std::array<u32,    N_PAIRS> _c0_x; // cluster peak strip index
	std::array<double, N_PAIRS> _e0_y; // accumulated cluster energy
	std::array<double, N_PAIRS> _delta_y;
	std::array<u32,    N_PAIRS> _c0_y; // cluster peak strip index
#endif

	RNFOOTTrack() = default;
	RNFOOTTrack(const std::array<double, 2>& , const std::array<double, 2>& , double , double, std::size_t
	#ifdef MND_FOOTTRACK_DEBUG
		,
		const std::array<double, N_PAIRS>& , // u
		const std::array<double, N_PAIRS>& , // v
		const std::array<double, N_PAIRS>& , // zx
		const std::array<double, N_PAIRS>& , // zy
		const std::array<double, N_PAIRS>& , // qx
		const std::array<double, N_PAIRS>& , // qy
		const std::array<double, N_PAIRS>& , // e0_x
		const std::array<double, N_PAIRS>& , // delta_x
		const std::array<u32,    N_PAIRS>& , // c0_x
		const std::array<double, N_PAIRS>& , // e0_y
		const std::array<double, N_PAIRS>& , // delta_y
		const std::array<u32,    N_PAIRS>&   // c0_y
	#endif
	);
	virtual ~RNFOOTTrack() = default;
	ClassDef(RNFOOTTrack, 1);
};
#ifndef MND_FOOTTRACK_DEBUG
static_assert(sizeof(RNFOOTTrack) == 64, "Alignment obsession, friends.");
#endif

struct RNFOOTHit {
	static constexpr u32 N_PAIRS = N_FOOT_DETECTORS / 2;
	struct Vertex {
		double x = NAN, y = NAN, z = NAN;

		Vertex() = default;
		Vertex(const mnd::geom::Point3D& p) : x(p.x), y(p.y), z(p.z) {}

		inline void Clean() noexcept { *this = Vertex{}; }

		virtual ~Vertex() = default;
		ClassDef(Vertex, 1);
	};

	std::array<RNFOOTPair, N_PAIRS> pair;
	Vertex vertex;
	std::vector<RNFOOTTrack> track;

    /* Count number of FOOT detector pairs with present both (X,Y) measurements. */
    u32 DataPresentCount() const noexcept;
    /* Returns true if all underlying FOOT pairs have a non-empty entry. */
    bool HasData() const noexcept;

	inline void Clean() noexcept {
		for(auto& p: pair) { p.Clean(); }
		vertex.Clean();
		track.clear();
	}

	RNFOOTHit() = default;
	virtual ~RNFOOTHit() = default;
	ClassDef(RNFOOTHit, 1);
};

template<size_t> struct TrackingGeometry;

struct SFOOTLayerGeometry {
	std::array<double, 8> A; // explicitly flattened row-major
	std::array<double, 2> dxy;
	std::array<double, 2> z;
	SFOOTLayerGeometry() = default;

	virtual ~SFOOTLayerGeometry() = default;
	ClassDef(SFOOTLayerGeometry, 1);
};

struct STrackingGeometry {
	std::array<SFOOTLayerGeometry, RNFOOTHit::N_PAIRS> layer;
	STrackingGeometry() = default;

	virtual ~STrackingGeometry() = default;
	ClassDef(STrackingGeometry, 1);
};

std::unique_ptr<TrackingGeometry<RNFOOTHit::N_PAIRS>> tg_from_stg(const STrackingGeometry& );

struct TFOOTHitCont : TContainer<RNFOOTHit> {
	static constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;

	/* Convention. In Kalman coordinates, place target's FOOT-facing edge nominally at 0.0.
	 * Will be shifted back to "real" FRS coordinates later. */
	static constexpr double TARGET_Z = 0.0;
	
	FOOTBoxParam* box;
	std::array<FOOTParam, 2>* foot_param[N_PAIRS];
	STrackingGeometry* geometry;

	std::string* setupName;

	std::array<double,3>* cost_coeff;
	TParameter<double> *max_cost, *max_cost_f;

	TH1I* h1_qtrack;
	TH1I* h1_track_nsampled;

	TH1I* diff_heavy_frag_vs_upstream;

#ifdef MND_FOOTTRACK_DEBUG
	/* Few histograms to allow for logging subsequent "innovations" in the Kalman procedure. */
	TH1I *h1_diff_q[N_PAIRS]; // Innovation in the difference: Qpair - Qpredicted
	TH1I *h1_diff_r[N_PAIRS]; // Innovation in the difference: |vec(r)_pair - vec(r)_predicted|
	TH1I *h1_diff_t[N_PAIRS]; // Innovation in the projected track's matching to upstream at target (distance)
	
	TH2I *h2_diff_q[N_PAIRS]; // Innovation in the difference: Qpair - Qpredicted *versus* (Qpair(x) - Qpair(y))
	TH2I *h2_diff_r[N_PAIRS]; // Innovation in the difference: vec(r)_pair - vec(r)_predicted
	TH2I *h2_diff_t[N_PAIRS]; // Innovation in the projected track's matching to upstream at target (vector)
#endif

	TFOOTHitCont();

	void Setup() override;
	void Init(TDictInfo ) override;
};

extern template struct HitMatrix<RNFOOTPair>; // instantiated in TFOOTHitProc.cxx
extern template struct Track<TFOOTHitCont::N_PAIRS + 1, RNFOOTPair>; // instantiated in TFOOTHitProc.cxx

/* Convert to 3D line. Delegates to the underlying deref operator. */
mnd::geom::Line3D RNTrackToLine3D(const RNFOOTTrack& ) noexcept;
