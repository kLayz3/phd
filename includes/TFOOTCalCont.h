#pragma once

#include "TContainer.hxx"
#include "TFOOTPedestalCont.h"
#include <cstddef>
#include <tuple>
#include <type_traits>

class TH1D;

struct TFOOTCluster {
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
	auto&& get() & { return get_helper<I>(*this); }

	template<std::size_t I>
	auto&& get() && { return get_helper<I>(*this); }

	template<std::size_t I>
	auto&& get() const & { return get_helper<I>(*this); }

	template<std::size_t I>
	auto&& get() const && { return get_helper<I>(*this); }

	TFOOTCluster(double x, double e, double m, ClusterType t) :
		fCX(x), fCE(e), fCM(m), fCT(t) {}
	virtual ~TFOOTCluster() = default;
	ClassDef(TFOOTCluster, 1);

private:
	template<std::size_t I, typename T>
	auto&& get_helper(T&& t) {
		static_assert(I<3, "Index out of bounds for TFOOTCluster");
		if constexpr(I == 0) 
			return std::forward<T>(t).fCX;
		else if constexpr(I == 1) 
			return std::forward<T>(t).fCE;
		else 
			return std::forward<T>(t).fCM;
	}
};

/* Hacks to make it structured-binding decomposable. */
namespace std {
	template<> struct tuple_size<::TFOOTCluster> : integral_constant<size_t, 3> {};
	template<> struct tuple_element<0, ::TFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, ::TFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, ::TFOOTCluster> { using type = double; };
}

struct RNFOOTCalCont {
	static constexpr size_t INIT_CAPACITY = _FOOT_N_STRIPS_PER_ASIC;
	static constexpr int N_STRIPS = _FOOT_N_STRIPS;

	using ClusterType = TFOOTCluster::ClusterType;
	enum Orientation {
		kUNSPECIFIED = 0,
		kX           = 1,
		kY           = 2,
	};

	Orientation o{};
	std::vector<TFOOTCluster> fCl{}; 

	/* Record whole event in a vector, if we find a large cluster, for some reason. */
	std::vector<double> _fBadE{};      /* Size will be either 0 or 640. */
	std::vector<double> _fHeClSize1{}; /* Size will be either 0 or 640. */
	
	RNFOOTCalCont() { 
		fCl.reserve(INIT_CAPACITY); 
		_fBadE.reserve(N_STRIPS);
		_fHeClSize1.reserve(N_STRIPS);
	}
	
	inline void Clean() noexcept { 
		fCl.clear(); 
		_fBadE.clear();
		_fHeClSize1.clear();
	}
	inline std::vector<double> E() const noexcept {
		std::vector<double> res;
		res.reserve(fCl.size());
		for(auto const& c : fCl) res.push_back(c.fCE);
		return res;
	}
	inline std::vector<double> X() const noexcept {
		std::vector<double> res;
		res.reserve(fCl.size());
		for(auto const& c : fCl) res.push_back(c.fCX);
		return res;
	}
	void AddCluster(TFOOTCluster cl) {
		fCl.push_back(std::move(cl));
	}
	void AddCluster(double x, double e, double m, ClusterType ty) {
		fCl.emplace_back(x, e, m, ty);
	}
	virtual ~RNFOOTCalCont() = default;
	ClassDef(RNFOOTCalCont, 1);
};

class TFOOTCalCont : public TContainer<RNFOOTCalCont> {
	friend class TFOOTCalProc;
public:	
	static constexpr int N_STRIPS = RNFOOTCalCont::N_STRIPS;
	using Orientation = RNFOOTCalCont::Orientation;	
	
	Orientation o = RNFOOTCalCont::kUNSPECIFIED; 
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

	TFOOTCalCont();

	void Init(TDictInfo info) /* override */;
};
