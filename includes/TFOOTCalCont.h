#pragma once

#include "monad/monad.hxx"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>

#include "TFOOTMapCont.h"

#include "util/json_struct_def.hh"
#include "util/PolyFitter.h"
#include "util/Interpolator.hxx"

#include "ExperimentAssumptions.hh"

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"

class TH1D;

enum class Orientation { X, Y, UNKNOWN };

struct FMultiPoly {
	using Vec = std::vector<double>;
	
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(i32, Z,   -1, 0);
	ADD_SERIALIZABLE_FIELD(Vec, pol, {}, 1);

	inline bool operator<(const FMultiPoly& rhs) const noexcept { return Z < rhs.Z; }
	FMultiPoly() = default;
	FMultiPoly(i32 Z_) : Z(Z_), pol{} {};

	virtual ~FMultiPoly() = default;
	ClassDef(FMultiPoly, 1);
};
ADD_JSON_TYPE_RESOLUTION(FMultiPoly, 1);

struct FOOTAsicGainParam {
	using Vec = std::vector<FMultiPoly>;
	GET_HELP_AUX_IMPL;
	
	/* Ordered vector, in `Z`. If not, analysis is fucked. Soothe paranoia w/ `IsSane()` method. */
	ADD_SERIALIZABLE_FIELD(Vec, multi_poly, {}, 0);

	const FMultiPoly* GetPoly(int Z) const noexcept;
	FMultiPoly* GetPoly(int Z) noexcept;

	/* Get the correct polynomial for Z, and evaluate it with x. If polynomial not found, is NAN. */
	double GetReferentMeasurement(int Z, double x) const;

	/* Check that poly's are sorted according to the Z attribute. */
	bool IsSaneZ() const;

	std::vector<i32> GetZ() const;

	/* Check for each 'x' that the each consecutive nominal energy evaluation along the `multi_poly` yields decreasing result.
	 * AKA: that for some 'x', the nominal value of Z cannot be above the value of Z-1. */
	bool IsSane(double x, std::vector<double>& evals) const noexcept;

	/* Return an iterator to the first referent measurment line (along E; for set x) that has ref >= e */
	auto GetMeasurementBound(double x, double e) const -> typename Vec::const_iterator;

	/* Return number of measurements, aka the size of the underlying multi_poly vector. */
	inline auto GetNReferentMeasurements() const noexcept { return multi_poly.size(); }
	
	FOOTAsicGainParam() = default;
	virtual ~FOOTAsicGainParam() = default;
	ClassDef(FOOTAsicGainParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTAsicGainParam, 0);

struct FOOTReferentADCMeasurement {
	GET_HELP_AUX_IMPL;
	ADD_SERIALIZABLE_FIELD(i32,    Z,   0,   0);
	ADD_SERIALIZABLE_FIELD(double, val, NAN, 1);

	inline bool operator<(const FOOTReferentADCMeasurement& rhs) const noexcept {
		return (Z < rhs.Z) and (val < rhs.val);
	}
	inline bool Valid() const noexcept {
		return Z > 0 && Z <= (i32)mnd::assume::primary.Z
			&& std::isfinite(val) && val > 0.0;
	}
	
	FOOTReferentADCMeasurement() = default;
	virtual ~FOOTReferentADCMeasurement() = default;
	ClassDef(FOOTReferentADCMeasurement, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTReferentADCMeasurement, 1);

struct FOOTGainParam {
	static constexpr int N_STRIPS_PER_ASIC = TFOOTMapCont::N_STRIPS_PER_ASIC;
	using AsicArray = std::array<FOOTAsicGainParam, TFOOTMapCont::N_ASIC>;
	using NominalValues = std::vector<FOOTReferentADCMeasurement>;

	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(AsicArray,     fit,           {}, 0);
	ADD_SERIALIZABLE_FIELD(NominalValues, nominal_value, {}, 1);

	/* Throws on invalid access. */
	inline const FOOTAsicGainParam& GetASIC(double x) const { return fit.at( x / N_STRIPS_PER_ASIC ); }

	std::vector<i32> GetNominalZ() const noexcept;
	
	MND_NEVER_INLINE void IsSaneThrowable() const;
	bool IsSane() const noexcept;

	void Reset() const noexcept;

	double GetNominalValue(i32 Z) const noexcept;

	/* Main call: get the correction factor g for the accumulated cluster energy e,
	 * based on the position x. Then the corrected energy is: e' = g*e */
	double CorrectionFactor(double x, double e) const;

	double Q(double E) const noexcept;

	/* Small helper function to plot what we're actually gain matching upon.
	 * Highly inefficient, but it's w/e. */
	[[ nodiscard ]] TH2D* GetHisto (
		int nbins_x = 640, int nbins_y = 500,
		int lo_y = 0, int hi_y = 2500
	) const;

	/* Get the 1D graph for a specific x-position. */
	[[ nodiscard ]] std::pair<TGraph*, TGraph*> GetGraph (
		const double x, int nbins_y = 500,
		int lo_y = 0, int hi_y = 4000
	) const;

	[[ nodiscard ]] std::pair<TGraph*, TLine*> GetRefZGraph(
		int Z, int Npts = 640
	) const;

	[[ nodiscard ]] std::vector<std::tuple<i32, TGraph*, TLine*>>
	GetAllRefZGraph(int Npts = 640) const;

	/* Get the 1D graph for charge converter. Never returns invalid pointers, even
	 * for empty calibration. Can throw on initializing the fitter. */
	[[ nodiscard ]] std::pair<TGraph*, TGraph*> GetQGraph (
		int npts = 200
	) const;

	FOOTGainParam() = default;
	virtual ~FOOTGainParam() = default;
	ClassDef(FOOTGainParam, 1);

private:
	/* Some cached values for quick Q- calculation.
	 * NB: if the object is re-evaluted, the values *need* to be recomputed, but the default
	 * JSON propagator cannot know this. Meaning that `Reset()` has to be called manually. */
	using InterpolatedCharge = mnd::Fitter<f64, mnd::PowerLaw<f64>, mnd::PowerLaw<f64>>;
	using ChargeConverter = std::variant<
		mnd::InvalidFitter<f64>, // [0] measurements (Z,E) available.
		mnd::PowerLaw<f64>, // [1,2,3] measurements (Z,E). In case a single measurement, fix the exponent 0.5
		InterpolatedCharge // >= 4  measurements (Z,E).
	>;

	/* This is to convert from E to charge., e.g.: Q(E) power law. */
	struct fit_coeff {
		mnd::PowerLaw<double> raw_fit;
		ChargeConverter converter;

		// add default ctor so it leaves the other field uninitialized
		fit_coeff() : raw_fit{NAN, NAN} {}
	};
	mutable mnd::cache<fit_coeff> fit_; //!
	mutable mnd::mutex mtx_; //!

	void Init() const;
};
ADD_JSON_TYPE_RESOLUTION(FOOTGainParam, 1);

struct FOOTDeltaFFT {
	using Coeff = std::array<double, 2>;
	using Vec = std::vector<Coeff>;
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(int, n, 0,  0);
	ADD_SERIALIZABLE_FIELD(Vec, c, {}, 1);

	double Evaluate(const double delta) const;
	
	/* Small helper function to plot what we're actually matching upon. */
	[[ nodiscard ]] TGraph* GetGraph(int Npts = 60) const;

	FOOTDeltaFFT() = default;
	virtual ~FOOTDeltaFFT() = default;
	ClassDef(FOOTDeltaFFT, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTDeltaFFT, 1);

/* Small struct to hold delta-correction associated params. */
struct FOOTDeltaParam {
	GET_HELP_AUX_IMPL;
	ADD_SERIALIZABLE_FIELD(double,       s,  0.33, 0);
	ADD_SERIALIZABLE_FIELD(double,       f,  0.50, 1);
	ADD_SERIALIZABLE_FIELD(FOOTDeltaFFT, f2, {},   2);

	/* Factor in range <0,1] corresponding to the initial triangular spline. */
	double CorrectionBasic(double delta) const noexcept;

	/* Total factor. Corrected cluster energy is then: e' = e / CorrectionFactor(δ) */
	double CorrectionFactor(double delta) const noexcept;
	
	/* Small helper function to plot what we're actually matching upon. */
	[[ nodiscard ]] TGraph* GetGraph(int Npts = 60, double x_lo = -0.5, double x_hi = 0.5) const;

	/* Small helper function to plot what we're actually matching upon.
	 * First argument is the "central ADC" gain-matched value. */
	[[ nodiscard ]] TGraph* GetGraphBasic(double E0 = 1.0, int Npts = 60) const;

	FOOTDeltaParam() = default;
	virtual ~FOOTDeltaParam() = default;
	ClassDef(FOOTDeltaParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTDeltaParam, 2);

struct RNFOOTCluster;

/* Parameters specific to a single FOOT detector. */
struct FOOTParam {
	constexpr static int N_STRIPS = TFOOTMapCont::N_STRIPS; // 640
	constexpr static double DETECTOR_MIDPOINT = static_cast<double>(N_STRIPS - 1) / 2; // 319.5 strip units
	constexpr static double DETECTOR_SIZE     = 96.0; // [ mm ]
	constexpr static double STRIP_TO_MM       = DETECTOR_SIZE / N_STRIPS; // 0.150 mm/strip
	constexpr static double MM_TO_STRIP       = 1.0 / STRIP_TO_MM; // 6.67 strip/mm

	constexpr static double CENTRE_THR_DEFAULT = 3.3;
	constexpr static double NEIGHB_THR_DEFAULT = 1.3;

	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(i32,              N,           -1,                 0);
	ADD_SERIALIZABLE_FIELD(double,           dz,          0.0,                1);
	ADD_SERIALIZABLE_FIELD(std::string,      orientation, {},                 2);
	ADD_SERIALIZABLE_FIELD(i32,              mirrored,    -1,                 3);
	ADD_SERIALIZABLE_FIELD(double,           c_threshold, CENTRE_THR_DEFAULT, 4);
	ADD_SERIALIZABLE_FIELD(double,           n_threshold, NEIGHB_THR_DEFAULT, 5);
	ADD_SERIALIZABLE_FIELD(double,           delta_p,     0.0,                6);
	ADD_SERIALIZABLE_FIELD(double,           delta_a,     0.0,                7);
	ADD_SERIALIZABLE_FIELD(FOOTGainParam,    gain,        {},                 8);
	ADD_SERIALIZABLE_FIELD(FOOTDeltaParam,   de,          {},                 9);

	Orientation GetOrientation() const noexcept;
	double R() const noexcept;

	/* Convert raw cluster position `cx` and its integrated ADC value
	 * `ce` into the final calibrated ADC value. */
	double E(const RNFOOTCluster& ) const noexcept;

	/* Main method: convert cluster energy (E) to nominal charge (Q)
	 * This energy E, has to be properly both gain-matched and delta-corrected. */
	double Q(double E) const noexcept;
	double Q(const RNFOOTCluster& ) const noexcept;

	/* Calculate the (bare) cluster hit position in millimeters.
     * Note: this does NOT take offsets or rotations into account. */
	double BarePosition(const RNFOOTCluster& ) const noexcept;

	/* Calculate the cluster hit position in millimeters, along the detector's axis. */
	double X0(const RNFOOTCluster& ) const noexcept;

	/* Convert bare measurement to strip number. Inverse of BarePosition. */
	double BareToStripnum(double ) const noexcept;

	/* Convert already offsetted measurement to strip number. Inverse of X0. */
	double X0ToStripnum(double ) const noexcept;

public:
	int de10_index_ = -1;
	virtual ~FOOTParam() = default;
	ClassDef(FOOTParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTParam, 9)

struct ExpertTarget {
	GET_HELP_AUX_IMPL
	ADD_SERIALIZABLE_FIELD(double, thickness, 0.0, 0)
	ADD_SERIALIZABLE_FIELD(double, width_x,   0.0, 1)
	ADD_SERIALIZABLE_FIELD(double, width_y,   0.0, 2)
	ADD_SERIALIZABLE_FIELD(double, dx,        0.0, 3)
	ADD_SERIALIZABLE_FIELD(double, dy,        0.0, 4)
	ADD_SERIALIZABLE_FIELD(double, dz,        0.0, 5)

	/* Target width in units of [mm] */
	double Width() const noexcept {
		return 10.0 * thickness / mnd::assume::s2::target_density;
	}

	virtual ~ExpertTarget() = default;
	ClassDef(ExpertTarget, 1);
};
ADD_JSON_TYPE_RESOLUTION(ExpertTarget, 5)

/* Parameters describing the whole FOOT box. Whatever the box may be :) */
struct FOOTBoxParam {
	GET_HELP_AUX_IMPL
	using Arr1 = std::array<double, N_FOOT_DETECTORS>;

	ADD_SERIALIZABLE_FIELD(double, z0,          NAN, 0)
	ADD_SERIALIZABLE_FIELD(double, width_outer, NAN, 1)
	ADD_SERIALIZABLE_FIELD(Arr1,   det_pos,     {},  2)
	ADD_SERIALIZABLE_FIELD(double, dx,          NAN, 3)
	ADD_SERIALIZABLE_FIELD(double, dy,          NAN, 4)
	ADD_SERIALIZABLE_FIELD(double, da,          NAN, 5)
	ADD_SERIALIZABLE_FIELD(double, db,          NAN, 6)
	ADD_SERIALIZABLE_FIELD(ExpertTarget, target, {}, 7)

	inline double GetTargetZ() const noexcept { return z0 - (width_outer / 2) + target.dz; }
	
	inline std::array<std::array<double, 2>, 2> TargetXYGeom() const noexcept {
		return {
			std::array {
				this->dx + target.dx - target.width_x / 2,
				this->dx + target.dx + target.width_x / 2,
			},
			std::array {
				this->dy + target.dy - target.width_y/ 2,
				this->dy + target.dy + target.width_y/ 2,
		}
		};
	}
	/* Calculate the absolute z- position of n-th FOOT detector in the box. */
	inline double GetFOOTZ(const int n, const FOOTParam* p = nullptr) const noexcept {
		double zfoot = z0 - (width_outer / 2) + det_pos.at(n);
		if(p) {
			zfoot += p->dz;
		}
		return zfoot;
	}
	inline double GetFOOTZ(const FOOTParam* p) const noexcept { return GetFOOTZ(p->N, p); }

	/* Calculate z- position of n-th FOOT detector in the box, relative to the target */
	inline double GetFOOTZRel(const int n, const FOOTParam* p = nullptr) const noexcept {
		return GetFOOTZ(n, p) - GetTargetZ();
	}
	inline double GetFOOTZRel(const FOOTParam* p) const noexcept {
		return GetFOOTZ(p) - GetTargetZ();
	}
	
	virtual ~FOOTBoxParam() = default;
	ClassDef(FOOTBoxParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTBoxParam, 7)

/* f(x; (a0,mu,sigma)) = a0 * exp( -0.5 * ((x-mu)/sigma)^2 ) */
struct FOOTClusterFit {
	static constexpr double TWO_PI = 6.283185307179586;
	static constexpr double TWO_PI_SQRT = 2.5066282746310002;
	static constexpr double TWO_PI_SQUARED = 19.739208802178716;
	static constexpr u32 LARGE_CLUSTER_CUTOFF = 7;

	inline static double ffourier(u32 k, double sigma, double delta) {
		return exp(-TWO_PI_SQUARED * k*k * sigma*sigma) * cos(TWO_PI*k*delta);
	}

	double a0     = NAN; /* Fitted Gauss amplitude. */
	double mu     = NAN; /* Fitted Gauss mean. */
	double sigma  = NAN; /* Fitted gauss width. */
	double delta  = NAN; /* Offset between mean and max strip. Also a null indicator. */
	
	inline double X() const noexcept { return mu; };
	inline int   I0() const noexcept { return mnd::rround<int>( mu - delta ); }
	inline double E() const noexcept { return a0 * sigma * TWO_PI_SQRT; }
	
	inline double E_discrete() const noexcept {
		return this->E() * (1 + 2 * (ffourier(1,sigma,delta) + ffourier(2,sigma,delta) + ffourier(3,sigma,delta)));
	}

	/* If true, all fit parameters are given. */
	bool IsOk() const noexcept {
		return mnd::isfinite(a0, mu, sigma, delta) &&
			sigma > 0 && a0 > 0;
	}

	FOOTClusterFit() = default;
	virtual ~FOOTClusterFit() = default;
	ClassDef(FOOTClusterFit, 1);
};

struct RNFOOTCluster {
	GET_HELP_AUX_IMPL
	enum ClusterType : u32 {
		kUNKNOWN    = 0, /* Unqualified. */
		kGOOD       = 1, /* Good cluster. Monotonically rising ADC values to peak ADC strip, then monotonically decreasing. */
		kFRAGMENTED = 2, /* One strip in noise; is missing between two sequences of the cluster. */
		kWAVEY      = 3, /* Cluster such as: _/\/\_, initial sequence isn't strictly rising, latter isn't stricly falling. */
		kMERGED     = 4, /* When two or more non-neighbouring distinct strips pass C-threshold check, and form a cluster. */
		kISOLATED   = 5, /* When only a single strip passes C-threshold */
		kHEAVY_ISOLATED = 6 /* When only a single strip passes C-threshold, far above the noise, Say 10σ above the noise. */
	};

	double fCX = 0; /* Cluster mean strip position. */
	double fCE = 0; /* Cluster summed energy. */
	u32    fCM = 0; /* Cluster multiplicity. */
	ClusterType fCT{}; /* Cluster type. */

	inline double Delta() const noexcept { return mnd::rround<int>(fCX) - fCX; }

	RNFOOTCluster(double, double, u32, ClusterType);
	RNFOOTCluster() = default;
	virtual ~RNFOOTCluster() = default;
	ClassDef(RNFOOTCluster, 1);

private:
	template<std::size_t I, typename Self>
	static decltype(auto) get_helper(Self&& self) noexcept {
		if constexpr(I == 0)      return (std::forward<Self>(self).fCX);
		else if constexpr(I == 1) return (std::forward<Self>(self).fCE);
		else if constexpr(I == 2) return (std::forward<Self>(self).fCM);
		else if constexpr(I == 3) return (std::forward<Self>(self).fCT);
		else static_assert(I < 4, "Index out of bounds for RNFOOTCluster::get");
	}
};

/* Make it structured-binding decomposable. */
namespace std {
	template<> struct tuple_size<RNFOOTCluster> : integral_constant<size_t, 4> {};
	template<> struct tuple_element<0, RNFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, RNFOOTCluster> { using type = u32; };
    template<> struct tuple_element<3, RNFOOTCluster> { using type = RNFOOTCluster::ClusterType; };
}

struct RNFOOTCal {
	static constexpr size_t INIT_CAPACITY = _FOOT_N_STRIPS_PER_ASIC;
	static constexpr int N_STRIPS = _FOOT_N_STRIPS;

	using ClusterType = RNFOOTCluster::ClusterType;
	std::vector<RNFOOTCluster> fCl{};

    /* Record whole event in a vector, if we find a large cluster, for some reason, or
     * a good cluster but with some weird gimmicks. Vector size will be either 0 or 640.*/
	std::vector<double> fRaw{};

    /* Type of currently captured event. */
    enum EventType : u32 {
        Undetermined,    // Not handled by the analysis yet
        Empty,           // Event doesn't contain FOOT data
        MassiveCluster,  // Event contains a massive cluster
        IsolatedStrip,   // Event contains an isolated high-adc strip (cluster size 1)
        Good,            // Otherwise, a good sanity-checked event
    } w_ev_type { EventType::Undetermined };

	RNFOOTCal();
	void Clean() noexcept;

	/* Return a fresh vector containing all the collected energies in the event. */
	std::vector<double> E() const noexcept;

	/* Return a fresh vector containing all the collected cluster positions in the event.
	 * It is in units of strips. */
	std::vector<double> X() const noexcept;
	inline void AddCluster(RNFOOTCluster cl) noexcept {
		fCl.push_back(std::move(cl));
	}
	inline void AddCluster(double x, double e, u32 m, ClusterType ty) noexcept {
		fCl.emplace_back(x, e, m, ty);
	}
    inline bool IsEmptyEvent() const noexcept {
        return w_ev_type == EventType::Empty;
    }

    /* A check if the debug field isn't already assigned. No need to assign it multiple times. */
    inline bool DebugAssigned() const noexcept {
        return !fRaw.empty();
    }
	virtual ~RNFOOTCal() = default;
	ClassDef(RNFOOTCal, 1);
};
static_assert(sizeof(RNFOOTCal) == 64);

inline void Add(FOOTParam&, const FOOTParam&) {}
inline void Add(FOOTBoxParam&, const FOOTBoxParam&) {}

struct TFOOTCalCont  : TContainer<RNFOOTCal> {
	friend struct TFOOTCalProc;

	constexpr static int N_STRIPS = RNFOOTCal::N_STRIPS;
	constexpr static int N_ASIC = TFOOTMapCont::N_ASIC;
	
	int FOOT_N = -1; /* Comes from sort step, isn't in order. */

	TH1I* h1_mult;
	TH1I* h1_dE;
	TH1I* h1_X;
	TH1I* h1_cl_type;

	TH2I* h1_dE_m3_light;
	TH2I* h1_dE_m3_heavy;

	TH2I* h2_mult_e;

	FOOTParam* setup;
	FOOTBoxParam* box;      // Optional; not every FOOT will register the box/name objects
	std::string* setupName; // based on `should_register_box_` predicate.

	TFOOTCalCont() = default;

	void Init(TDictInfo info) override;
	void Setup() override;
	inline void SetRegisterBox(bool b) { should_register_box_ = b; }

private:
	FOOTParam par;   /* Local object, will just get copied around. Is fine. */
	FOOTBoxParam bpar; /* Again, local object for a small temporary buffer. */
	bool should_register_box_ = false;
	
	/* ^^^ I don't make them internally linked in the .cxx since different containers
	 * could point to different parameter files/objects,.. then I'd need some map and query it..
	 * It's a bit too complex for less than 500 byte of memory saved. */

};
