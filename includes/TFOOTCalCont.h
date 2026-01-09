#pragma once

#include "monad/monad.hxx"
#include <cstddef>

#include "TFOOTMapCont.h"

#define GET_HELP_AUX_IMPL  \
	template<std::size_t I> \
    decltype(auto) get() &       noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() const & noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() &&      noexcept { return get_helper<I>(std::move(*this)); } \
    template<std::size_t I> \
    decltype(auto) get() const&& noexcept { return get_helper<I>(std::move(*this)); }

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
	decltype(auto) get() &        noexcept { return get_helper<I>(*this); }

	template<std::size_t I>
	decltype(auto) get() const &  noexcept { return get_helper<I>(*this); }

	template<std::size_t I>
	decltype(auto) get() &&       noexcept { return get_helper<I>(std::move(*this)); }

	template<std::size_t I>
	decltype(auto) get() const && noexcept { return get_helper<I>(std::move(*this)); }

	RNFOOTCluster(double x, double e, double m, ClusterType t);
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
	inline void AddCluster(double x, double e, double m, ClusterType ty) noexcept {
		fCl.emplace_back(x, e, m, ty);
	}
	virtual ~RNFOOTCal() = default;
	ClassDef(RNFOOTCal, 1);
};

struct TFOOTCalCont : TContainer<RNFOOTCal> {
	friend struct TFOOTCalProc;

	constexpr static int N_STRIPS = RNFOOTCal::N_STRIPS;
	constexpr static double CENTRE_THR_DEFAULT = 3.3;
	constexpr static double NEIGHB_THR_DEFAULT = 1.3;

	/* These values will not get serialized. */
	int FOOT_N = -1; /* Comes from sort step. */
	int POS = -1; /* Force the FOOT's to be labelled 0,1,2,3,4,5,6,7 from now onward. */

	double c_threshold = CENTRE_THR_DEFAULT;
	double n_threshold = NEIGHB_THR_DEFAULT;

	TH1I* h1_raw_mult; 
	TH1I* h1_mult; 
	TH1I* h1_dE; 
	TH1I* h1_X; 
	TH1I* h1_cl_type; 

	TH1I* h1_dE_m1; 
	TH1I* h1_dE_m2; 
	TH1I* h1_dE_m3; 
	TH1I* h1_sn_ratio; 

	std::array<double,2> *threshold;

	TFOOTCalCont() = default;

	void Init(TDictInfo info) override;
	void Setup() override;
};
