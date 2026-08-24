#pragma once

#include "monad/monad.hxx"
#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "TFOOTCalCont.h"

class TH2I;

/* Heavy debugging features compiled in. Can be toggled together automatically.
 * This is either for sanity debugging or for: ../scripts/hit/foot_track_analysis.C */
//#define MND_FOOTTRACK_DEBUG

/* In case the heavy debug build is enabled, also feature it in here. */
#if defined(MND_DEBUG_ENABLED) && !defined(MND_FOOTTRACK_DEBUG)
#	define MND_FOOTTRACK_DEBUG
#endif

#ifdef MND_FOOTTRACK_DEBUG
#	ifndef MND_HITMATRIX_DO_BOUNDS_CHECK
#		define MND_HITMATRIX_DO_BOUNDS_CHECK
#	endif
#endif

#include "util/HitMatrix.hxx"
#include "util/FTrack.h"

/* Represent the 'charge' measurement of each layer. */
struct FOOTQ {
	using ClusterType = RNFOOTCluster::ClusterType;

	float q; // nominal charge 'value'
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
	double m; // measurement [ mm ]
	
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

    /* Nominally, it is the z- position of the pair center in FOOT-box coordinates.
     * However, we can use the lowest two bits of the mantissa to indicate if the
     * `x` or the `y` foot are missing. This changes the precision by ~1e-16 or so
     * such that it doesn't really matter. Sane person would implement a separate flag, but I'm an
     * OCD masochist who likes alignment very much. Not safe for kids. :)
     * NOTE: Remapping it to `nan` would open up the full 51-bits mantissa, but then have to handle it properly per event. */
	double z;
	RNFOOTPair() = default;

    constexpr static u32 _DATA_X_PRESENT_BITINDEX = 0;
    constexpr static u32 _DATA_Y_PRESENT_BITINDEX = 1;
    constexpr static u64 _DATA_X_PRESENT_BITMASK  = (1ULL << _DATA_X_PRESENT_BITINDEX);
    constexpr static u64 _DATA_Y_PRESENT_BITMASK  = (1ULL << _DATA_Y_PRESENT_BITINDEX);

    bool HasDataX() const noexcept;
    bool HasDataY() const noexcept;
    bool HasData() const noexcept;

	inline void Clean() noexcept { x.clear(); y.clear(); }
	virtual ~RNFOOTPair() = default;
	ClassDef(RNFOOTPair, 1);
};
static_assert(sizeof(RNFOOTPair) == 64, "Alignment obsession, friends.");

struct RNFOOTTrack {
	double x0, y0;
	double ax, ay;
	double Q;

	double score;
	std::size_t n = 0; // number of collected pts, size_t anyway chosen as it will be aligned to 8-byte

	inline bool operator<(const RNFOOTTrack& rhs) const noexcept {
		return Q > rhs.Q; // descending charge.
	}
#ifdef MND_FOOTTRACK_DEBUG
	static constexpr u32 N_PAIRS = N_FOOT_DETECTORS / 2;
	
	std::array<double, N_PAIRS> _x;
	std::array<double, N_PAIRS> _y;
	std::array<double, N_PAIRS> _z;
	std::array<double, N_PAIRS> _q;
	std::array<double, N_PAIRS> _sq;

	std::array<double, N_PAIRS> _e0_x;
	std::array<double, N_PAIRS> _delta_x;
	std::array<u32,    N_PAIRS> _c0_x; // cluster peak value
	std::array<double, N_PAIRS> _e0_y;
	std::array<double, N_PAIRS> _delta_y;
	std::array<u32,    N_PAIRS> _c0_y; // cluster peak value
#endif

	RNFOOTTrack() = default;
	RNFOOTTrack(const std::array<double, 2>& , const std::array<double, 2>& , double , double, std::size_t
#ifdef MND_FOOTTRACK_DEBUG
		,
		const std::array<double, N_PAIRS>& , 
		const std::array<double, N_PAIRS>& , 
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<u32,    N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& , 
		const std::array<u32,    N_PAIRS>&
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

	RNFOOTTrack heavy_fragment; // Mostly a debug field. The same track will be found in the vector (usually).
	std::array<RNFOOTPair, N_PAIRS> pair;
	Vertex vertex;
	std::vector<RNFOOTTrack> track;

	RNFOOTHit() = default;

	inline void Clean() noexcept { 
		heavy_fragment.n = 0;
		for(auto& p: pair) { p.Clean(); }
		vertex.Clean();
		track.clear(); 
	}
	virtual ~RNFOOTHit() = default;
	ClassDef(RNFOOTHit, 1);
};

struct TFOOTHitCont : TContainer<RNFOOTHit> {
	static constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	
	FOOTBoxParam* box;
	std::array<FOOTParam, 2>* foot_param[N_PAIRS];
	std::string* setupName;

	std::array<double,3>* cost_coeff;
	TParameter<double> *max_cost, *max_cost_f;

	TH1I* h1_qtrack;
	TH1I* h1_track_nsampled;

	TH1I* diff_heavy_frag_vs_upstream;

#ifdef MND_FOOTTRACK_DEBUG
	TH1I *h1_diff_q[N_PAIRS];
	TH1I *h1_diff_r[N_PAIRS];
	TH1I *h1_diff_t[N_PAIRS];
	TH1I *h1_acc_q[N_PAIRS];
	TH1I *h1_acc_r[N_PAIRS];
	TH1I *h1_acc_t[N_PAIRS];
#endif

	TFOOTHitCont();

	void Setup() override;
	void Init(TDictInfo ) override;
};

extern template struct HitMatrix<RNFOOTPair>; // instantiated in TFOOTHitProc.cxx
extern template struct Track<TFOOTHitCont::N_PAIRS + 1, RNFOOTPair>; // instantiated in TFOOTHitProc.cxx

mnd::geom::Line3D RNTrackToLine3D(const RNFOOTTrack& );
