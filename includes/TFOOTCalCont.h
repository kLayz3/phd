#pragma once

#include "monad/monad.hxx"
#include <cstddef>

#include "TFOOTMapCont.h"
#include "json_struct_def.hh"
#include "TParameter.h"

class TH1D;

#define FOOT_ID_0 10
#define FOOT_ID_1 17
#define FOOT_ID_2 19
#define FOOT_ID_3 20
#define FOOT_ID_4 22
#define FOOT_ID_5 25
#define FOOT_ID_6 23
#define FOOT_ID_7 21

inline constexpr i32 static_detectors[] = {
	FOOT_ID_0, // Gets mapped to FOOT0
	FOOT_ID_1, // Gets mapped to FOOT1
	FOOT_ID_2, // Gets mapped to FOOT2
	FOOT_ID_3, // Gets mapped to FOOT3
	FOOT_ID_4, // Gets mapped to FOOT4
	FOOT_ID_5, // Gets mapped to FOOT5
	FOOT_ID_6, // Gets mapped to FOOT6
	FOOT_ID_7  // Gets mapped to FOOT7
};
inline constexpr i32 N_FOOT = mnd::len(static_detectors);
static_assert(N_FOOT == N_FOOT_DETECTORS);

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
	bool IsOk() const noexcept { return std::isfinite(delta); }

	FOOTClusterFit() = default;
	virtual ~FOOTClusterFit() = default;
	ClassDef(FOOTClusterFit, 1);
};

struct RNFOOTCluster {
	enum ClusterType {
		kUNKNOWN    = 0, /* Unqualified. */
		kGOOD       = 1, /* Good cluster. Monotonically rising ADC values to peak ADC strip, then monotonically decreasing. */
		kFRAGMENTED = 2, /* One strip in noise; is missing between two sequences of the cluster. */
		kWAVEY      = 3, /* Cluster such as: _/\/\_, initial sequence isn't strictly rising, latter isn't stricly falling. */
		kMERGED     = 4, /* When two or more non-neighbouring distinct strips pass C-threshold check, and form a cluster. */
	};

	double fCX = 0; /* Cluster mean strip position. */
	double fCE = 0; /* Cluster summed energy. */
	u32    fCM = 0; /* Cluster multiplicity. */
	ClusterType fCT{}; /* Cluster type. */
	FOOTClusterFit fit{};

	inline double Delta() const noexcept { return mnd::rround<int>(fCX) - fCX; }

	template<std::size_t I>
	decltype(auto) get() &        noexcept { return get_helper<I>(*this); }

	template<std::size_t I>
	decltype(auto) get() const &  noexcept { return get_helper<I>(*this); }

	template<std::size_t I>
	decltype(auto) get() &&       noexcept { return get_helper<I>(std::move(*this)); }

	template<std::size_t I>
	decltype(auto) get() const && noexcept { return get_helper<I>(std::move(*this)); }

	RNFOOTCluster(double, double, u32, ClusterType, FOOTClusterFit);
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
	template<> struct tuple_size<::RNFOOTCluster> : integral_constant<size_t, 4> {};
	template<> struct tuple_element<0, ::RNFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, ::RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, ::RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<3, ::RNFOOTCluster> { using type = RNFOOTCluster::ClusterType; };
}

struct alignas(mnd::CL) RNFOOTCal {
	static constexpr size_t INIT_CAPACITY = _FOOT_N_STRIPS_PER_ASIC;
	static constexpr int N_STRIPS = _FOOT_N_STRIPS;

	using ClusterType = RNFOOTCluster::ClusterType;
	std::vector<RNFOOTCluster> fCl{}; 

	/* Record whole event in a vector, if we find a large cluster, for some reason. */
	std::vector<double> _fBadE{};      /* Size will be either 0 or 640. */
	std::vector<double> _fHeClSize1{}; /* Size will be either 0 or 640. */
	
	RNFOOTCal();
	void Clean() noexcept;

	/**
	 * Return a fresh vector containing all the collected energies in the event.
	 */
	std::vector<double> E() const noexcept; 

	/**
	 * Return a fresh vector containing all the collected energies in the event.
	 */
	std::vector<double> X() const noexcept; 
	inline void AddCluster(RNFOOTCluster cl) noexcept {
		fCl.push_back(std::move(cl));
	}
	inline void AddCluster(double x, double e, u32 m, ClusterType ty, FOOTClusterFit f) noexcept {
		fCl.emplace_back(x, e, m, ty, f);
	}
	virtual ~RNFOOTCal() = default;
	ClassDef(RNFOOTCal, 1);
};

/* Keep this struct a mirror of the setup file JSON.
 * We don't serialize nlohman::json as a raw dump directly into ROOT, as parsing it
 * then can be annoying. Just solve all the parsing right here, and then have a unique type
 * that can be queried directly. */

struct FOOTParam {
	constexpr static double CENTRE_THR_DEFAULT = 3.3;
	constexpr static double NEIGHB_THR_DEFAULT = 1.3;

	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(i32,         N,           -1,                 0);
	ADD_SERIALIZABLE_FIELD(double,      z0,          NAN,                1);
	ADD_SERIALIZABLE_FIELD(std::string, orientation, {},                 2);
	ADD_SERIALIZABLE_FIELD(i32,         mirrored,    -1,                 3);
	ADD_SERIALIZABLE_FIELD(double,      c_threshold, CENTRE_THR_DEFAULT, 4);
	ADD_SERIALIZABLE_FIELD(double,      n_threshold, NEIGHB_THR_DEFAULT, 5);
	ADD_SERIALIZABLE_FIELD(double,      delta_a,     0.0,                6);
	ADD_SERIALIZABLE_FIELD(double,      delta_p,     0.0,                7);

	virtual ~FOOTParam() = default;
	ClassDef(FOOTParam, 1);
};
ADD_STD_TYPE_RESOLUTION_(FOOTParam, 7)

struct TFOOTCalCont  : TContainer<RNFOOTCal> {
	friend struct TFOOTCalProc;

	constexpr static int N_STRIPS = RNFOOTCal::N_STRIPS;
	/* These values will not get serialized. */
	int FOOT_N = -1; /* Comes from sort step, isn't in order. */
	FOOTParam par;   /* Local object, will just get copied around. Is fine. */

	TH1I* h1_mult; 
	TH1I* h1_dE; 
	TH1I* h1_X; 
	TH1I* h1_cl_type; 

	TH1I* h1_dE_m1; 
	TH1I* h1_dE_m2; 
	TH1I* h1_dE_m3; 
	TH1I* h1_cl_sigma;
	TH2I* h2_mult_e;

	FOOTParam* setup;

	TFOOTCalCont() = default;

	void Init(TDictInfo info) override;
	void Setup() override;
};
