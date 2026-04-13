#pragma once

#include "monad/monad.hxx"
#include <cstddef>

#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "util/PolyFitter.hxx"
#include "util/Tracking.hxx"

#include "TGraph.h"

class TH1D;

struct FOOTAsicGainParam {
	using Vec = std::vector<double>;
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(Vec, central, {}, 0);
	ADD_SERIALIZABLE_FIELD(Vec, lateral, {}, 1);

	inline double ValueCentral(const double x) const noexcept { return poly::Eval(x, central); }
	inline double ValueLateral(const double x) const noexcept { return poly::Eval(x, lateral); }

	FOOTAsicGainParam() = default;
	virtual ~FOOTAsicGainParam() = default;
	ClassDef(FOOTAsicGainParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTAsicGainParam, 1);

struct FOOTGainParam {
	static constexpr double CARBON_ADC = 3600;
	using AsicArray = std::array<FOOTAsicGainParam, TFOOTMapCont::N_ASIC>;
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(double,    lat_avg, 0.0, 0);
	ADD_SERIALIZABLE_FIELD(double,    mid_avg, 0.0, 1);
	ADD_SERIALIZABLE_FIELD(AsicArray, fit,     {},  2);
	
	inline std::array<double, 2> GetCentralLateralMeanE(const double x) const noexcept {
		int i = static_cast<int>( x / TFOOTMapCont::N_STRIPS_PER_ASIC );
		if(i >= (int)fit.size())
			ERROR("GetCentralLateralMeanE: Requested cluster at %.2f, index=%d > %d\n",
				x, i, (int)fit.size());
		
		const auto& p = fit[i]; 
		return { p.ValueCentral(x), p.ValueLateral(x) };
	}

	/**
	 * Get two values, corresponding to the: 
	 * 1) needed gain at the completely central
	 *   9C hit, which will be gain-matched into `CARBON_ADC`.
	 * 2) needed gain at the lateral (non-bonded strip) 9C hit, which will
	 * be gain matched into `CARBON_ADC * S/M `. 
	 */
	inline std::array <
		std::array<double, 2>, 2 
	> GetCentralLateralGainValue(const double x) const noexcept {
		auto [central_e, lateral_e ] = this->GetCentralLateralMeanE(x);	
		return {{
			/* pos            gain           */
			{central_e, CARBON_ADC / central_e}, 
			{lateral_e, (lat_avg/mid_avg) * CARBON_ADC / lateral_e} 
		}};
	}

	/* Factor to divide the measured ADC value to fit the equal C line. */
	inline double CorrectionFactor(const double x, const double e) const noexcept {
		auto values = this->GetCentralLateralGainValue(x);

		std::array<double, 2> lin_fit = GetLine(values);
		
		/* Point where d/de( g(e) * e) == 0; meaning that the gain matching will reorder two ADC values.
		 * At this point, the linear curve becomes invalid. */
		double gain_invalid_after = -lin_fit[0] / (2 * lin_fit[1]);
		gain_invalid_after = (gain_invalid_after > 0) ? gain_invalid_after : DBL_MAX;
		double gain_value_at_invalid_point = lin_fit[0] / 2;

		if(x > gain_invalid_after) 
			return e * gain_value_at_invalid_point;
		else 
			return poly::Eval(e, lin_fit);
	}

	/* Small helper function to plot what we're actually gain matching upon. */
	[[ nodiscard ]] inline std::pair<TGraph*, TGraph*> GetGraph() const {
		const int N =  (int)TFOOTMapCont::N_STRIPS;
		TGraph* g_cen = new TGraph(N);	
		TGraph* g_lat = new TGraph(N);	
		for(int i=0; i<N; ++i) {
			double x0 = i + 0.5;
			auto [y_cen, y_lat] = this->GetCentralLateralMeanE(x0);
			g_cen->SetPoint(i, x0, y_cen);
			g_lat->SetPoint(i, x0, y_lat);
		}
		g_cen->SetLineColor(kRed);
		g_cen->SetLineWidth(3);
		g_lat->SetLineColor(kMagenta + 1);
		g_lat->SetLineWidth(3);
		return {g_cen, g_lat};
	}

	FOOTGainParam() = default;
	virtual ~FOOTGainParam() = default;
	ClassDef(FOOTGainParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTGainParam, 2);

struct FOOTDeltaFFT {
	using Coeff = std::array<double, 2>;
	using Vec = std::vector<Coeff>;
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(int, n, 0,  0);
	ADD_SERIALIZABLE_FIELD(Vec, c, {}, 1);

	double Evaluate(const double delta) const;
		
	/* Small helper function to plot what we're actually matching upon. */
	[[ nodiscard ]] inline TGraph* GetGraph(int Npts = 60) const {
		constexpr double x_lo = -0.5;
		constexpr double x_hi =  0.5;

		TGraph* g = new TGraph(Npts);	
		for(int i=0; i<Npts; ++i) {
			double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts;
			double y0 = this->Evaluate(x0);
			g->SetPoint(i, x0, y0);
		}
		g->SetLineColor(kRed);
		g->SetLineWidth(3);
		return g;
	}

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
	inline double CorrectionBasic(const double delta) const noexcept {
		double x = std::min( std::abs(delta), s );
		return 1 - (1-f) / s * x;
	}
	/* Total factor. Corrected cluster energy is then: e' = e / CorrectionFactor(δ) */
	inline double CorrectionFactor(const double delta) const noexcept {
		
		double corr1 = this->CorrectionBasic(delta); // Interval [f, 1]
		double corr2 = this->f2.Evaluate(delta); // ~approx. [-0.3, 0.3]

		return corr1 * (1 + corr2);
	}
	
	/* Small helper function to plot what we're actually matching upon. */
	[[ nodiscard ]] inline TGraph* GetGraph(int Npts = 60, double x_lo = -0.5, double x_hi = 0.5) const {
		TGraph* g = new TGraph(Npts);	
		for(int i=0; i<Npts; ++i) {
			double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts;
			double y0 = this->CorrectionFactor(x0);
			g->SetPoint(i, x0, y0);
		}
		g->SetLineColor(kRed);
		g->SetLineWidth(3);
		return g;
	}


	FOOTDeltaParam() = default;
	virtual ~FOOTDeltaParam() = default;
	ClassDef(FOOTDeltaParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTDeltaParam, 2);

/* Parameters specific to a single FOOT detector. */
struct FOOTParam {
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

	int de10_index_ = -1;
	virtual ~FOOTParam() = default;
	ClassDef(FOOTParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTParam, 9)

/* Parameters describing the whole FOOT box. Whatever the box may be :) */
struct FOOTBoxParam {
	GET_HELP_AUX_IMPL
	using Arr1 = std::array<double, N_FOOT_DETECTORS>;

	ADD_SERIALIZABLE_FIELD(double, z0,          NAN, 0)
	ADD_SERIALIZABLE_FIELD(double, width_inner, NAN, 1)
	ADD_SERIALIZABLE_FIELD(Arr1,   det_pos,     {},  2)
	ADD_SERIALIZABLE_FIELD(double, width_outer, NAN, 3)
	ADD_SERIALIZABLE_FIELD(double, dx,          NAN, 4)
	ADD_SERIALIZABLE_FIELD(double, dy,          NAN, 5)
	ADD_SERIALIZABLE_FIELD(double, da,          NAN, 6)
	ADD_SERIALIZABLE_FIELD(double, db,          NAN, 7)

	inline double GetTargetZ() const noexcept { return z0 - (width_outer / 2); }

	/* Calculate the absolute z- position of n-th FOOT detector in the box. */
	inline double GetFOOTZ(const int n, const FOOTParam* p = nullptr) const noexcept {
		double zf = GetTargetZ() + det_pos.at(n);
		if(p) {
			/* For sanity check. The `N` field supplied there must be factor 2- up to `n` 
			 * since FOOTs by convention go in order. */
			if(n != (p->N / 2))
				WARN("Fetching FOOT z-position, but from box requested n=%d, "
					"supplied FOOT param handle reads .N=%d (position=%d)\n", n, p->N, p->N/2);
			zf += p->dz;
		}
		return zf;
	}
	inline double GetFOOTZ(const FOOTParam* p) const noexcept { return GetFOOTZ(p->N/2, p); }

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
	bool IsOk() const noexcept { return std::isfinite(delta); }

	FOOTClusterFit() = default;
	virtual ~FOOTClusterFit() = default;
	ClassDef(FOOTClusterFit, 1);
};

struct RNFOOTCluster {
	enum ClusterType : u32 {
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
	double fCP = 0; /* Cluster max energy in a single strip. */
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

	RNFOOTCluster(double, double, u32, ClusterType, double, FOOTClusterFit);
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
		else if constexpr(I == 4) return (std::forward<Self>(self).fCP);
		else static_assert(I < 5, "Index out of bounds for RNFOOTCluster::get");
	} 
};

/* Make it structured-binding decomposable. */
namespace std {
	template<> struct tuple_size<::RNFOOTCluster> : integral_constant<size_t, 5> {};
	template<> struct tuple_element<0, ::RNFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, ::RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, ::RNFOOTCluster> { using type = u32; };
    template<> struct tuple_element<3, ::RNFOOTCluster> { using type = RNFOOTCluster::ClusterType; };
    template<> struct tuple_element<4, ::RNFOOTCluster> { using type = double; };
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
	inline void AddCluster(double x, double e, u32 m, ClusterType ty, double p, FOOTClusterFit f) noexcept {
		fCl.emplace_back(x, e, m, ty, p, f);
	}
	virtual ~RNFOOTCal() = default;
	ClassDef(RNFOOTCal, 1);
};

struct TFOOTCalCont  : TContainer<RNFOOTCal> {
	friend struct TFOOTCalProc;

	constexpr static int N_STRIPS = RNFOOTCal::N_STRIPS;
	
	int FOOT_N = -1; /* Comes from sort step, isn't in order. */

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
	FOOTBoxParam* box; // Optional; not every FOOT will register the box object 
					   // based on `should_register_box_` predicate.
	TParameter<bool>* gain_matched;

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
