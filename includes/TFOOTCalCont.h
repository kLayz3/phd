#pragma once

#include "TContainer.hxx"
#include "TFOOTMapCont.h"
#include <cstddef>
#include <tuple>
#include <type_traits>

class TH1D;

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
	double fCM = 0; /* Cluster mean multiplicity. */
	ClusterType fCT{}; /* Cluster type. */

	template<std::size_t I>
	inline auto&& get() & { return get_helper<I>(*this); }

	template<std::size_t I>
	inline auto&& get() && { return get_helper<I>(*this); }

	template<std::size_t I>
	inline auto&& get() const & { return get_helper<I>(*this); }

	template<std::size_t I>
	inline auto&& get() const && { return get_helper<I>(*this); }

	RNFOOTCluster(double x, double e, double m, ClusterType t);
	RNFOOTCluster() = default;
	virtual ~RNFOOTCluster() = default;
	ClassDef(RNFOOTCluster, 1);

private:
	template<std::size_t I, typename T>
	auto&& get_helper(T&& t);
};

/* Hacks to make it structured-binding decomposable. */
namespace std {
	template<> struct tuple_size<::RNFOOTCluster> : integral_constant<size_t, 3> {};
	template<> struct tuple_element<0, ::RNFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, ::RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, ::RNFOOTCluster> { using type = double; };
}

struct RNFOOTCal {
	static constexpr size_t INIT_CAPACITY = _FOOT_N_STRIPS_PER_ASIC;
	static constexpr int N_STRIPS = _FOOT_N_STRIPS;

	using ClusterType = RNFOOTCluster::ClusterType;
	enum Orientation {
		kUNSPECIFIED = 0,
		kX           = 1,
		kY           = 2,
	};

	Orientation o{};
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
	inline void AddCluster(double x, double e, double m, ClusterType ty) noexcept {
		fCl.emplace_back(x, e, m, ty);
	}
	virtual ~RNFOOTCal() = default;
	ClassDef(RNFOOTCal, 1);
};

class TFOOTCalCont : public TContainer<RNFOOTCal> {
	friend class TFOOTCalProc;
public:	
	static constexpr int N_STRIPS = RNFOOTCal::N_STRIPS;
	using Orientation = RNFOOTCal::Orientation;	
	
	Orientation o = RNFOOTCal::kUNSPECIFIED; 
	int FOOT_N = -1; /* Comes from sort step. */
	int POS = -1; /* Force the FOOT's to be labelled 0,1,2,3,4,5,6,7 from now onward. */

	TH1I* h1_raw_mult; 
	TH1I* h1_mult; 
	TH1I* h1_dE; 
	TH1I* h1_X; 
	TH1I* h1_cl_type; 

	TH1I* h1_dE_m1; 
	TH1I* h1_dE_m2; 
	TH1I* h1_dE_m3; 
	TH1I* h1_sn_ratio; 

	TFOOTCalCont() = default;

	void Init(TDictInfo info) /* override */;
};
