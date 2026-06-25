#pragma once

#include "monad/monad.hxx"
#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "TFOOTCalCont.h"

class TH2I;

/* Debugging features compiled in. Can be toggled together automatically */
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

	float q; // nominal 'value'

	/* Few fields taken from RNFOOTCluster.. */
	u32 fCM = 0; /* Cluster multiplicity. */ 
	ClusterType fCT{}; /* Cluster type. */
	
	FOOTQ() = default;
	FOOTQ(float q_, u32 fCM_, ClusterType fCT_) :
		q(q_), fCM(fCM_), fCT(fCT_) {}

	virtual ~FOOTQ() = default;
	ClassDef(FOOTQ, 1);
};

struct FOOTHit {
	using ClusterType = FOOTQ::ClusterType;

	FOOTQ Q; 
	double m; // measurement [ mm ]
	
	FOOTHit() = default;
	FOOTHit(double q_, u32 fCM_, ClusterType fCT_, double m_) :
		Q(q_, fCM_, fCT_), m(m_) {}

	virtual ~FOOTHit() = default;
	ClassDef(FOOTHit, 1);
};
std::ostream& operator<<(std::ostream& , const FOOTHit& ) noexcept;

struct RNFOOTPair {
	std::vector<FOOTHit> x;
	std::vector<FOOTHit> y;
	double z;
	RNFOOTPair() = default;

	inline void Clean() noexcept { x.clear(); y.clear(); }
	virtual ~RNFOOTPair() = default;
	ClassDef(RNFOOTPair, 1);
};

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
#endif

	RNFOOTTrack() = default;
	RNFOOTTrack(const std::array<double, 2>& , const std::array<double, 2>& , double , double, std::size_t
#ifdef MND_FOOTTRACK_DEBUG
		,
		const std::array<double, N_PAIRS>& , 
		const std::array<double, N_PAIRS>& , 
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& ,
		const std::array<double, N_PAIRS>& 
#endif
	);
	virtual ~RNFOOTTrack() = default;
	ClassDef(RNFOOTTrack, 1);
};

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

	RNFOOTTrack heavy_fragment; // mostly a debug field. The same track will be found in the vector (usually).
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

	std::array<double,4>* cost_coeff;
	TParameter<double>* max_cost;

	TH1I* h1_qtrack;
	TH1I* h1_track_nsampled;

	TH1I* diff_heavy_frag_vs_upstream;

	TFOOTHitCont();

	void Setup() override;
	void Init(TDictInfo ) override;
};

extern template struct HitMatrix<RNFOOTPair>; // instantiated in TFOOTHitProc.cxx
extern template struct Track<TFOOTHitCont::N_PAIRS + 1, RNFOOTPair>; // instantiated in TFOOTHitProc.cxx

mnd::geom::Line3D RNTrackToLine3D(const RNFOOTTrack& );
