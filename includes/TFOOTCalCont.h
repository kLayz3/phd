#pragma once

#include "monad/monad.hxx"
#include <algorithm>
#include <cstddef>

#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "util/PolyFitter.hxx"
#include "util/Tracking.hxx"

#include "TGraph.h"

class TH1D;

using ExtrapolateLowZ = mnd::BinaryOpt;

struct FMultiPoly {
	using Vec = std::vector<double>;
	
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(i32, Z,   -1, 0);
	ADD_SERIALIZABLE_FIELD(Vec, pol, {}, 1);

	bool operator<(const FMultiPoly& rhs) const noexcept { return Z < rhs.Z; }
	FMultiPoly() = default;
	FMultiPoly(i32 Z_) : Z(Z_), pol{} {};

	virtual ~FMultiPoly() = default;
	ClassDef(FMultiPoly, 1);
};
ADD_JSON_TYPE_RESOLUTION(FMultiPoly, 1);

struct FOOTAsicGainParam {
	using Vec = std::vector<FMultiPoly>;
	GET_HELP_AUX_IMPL;
	
	static constexpr double PROTON_ADC = 100;
	inline static i32 NominalE(i32 Z) noexcept { return Z*Z * PROTON_ADC; }

	/* Ordered vector, in `Z`. If not, analysis is fucked. If paranoid, check w/ `IsSane()` method?*/
	ADD_SERIALIZABLE_FIELD(Vec, multi_poly, {}, 0);

	inline const FMultiPoly* GetPoly(int Z) const noexcept {
		for(const auto& p : multi_poly) { if(p.Z == Z) return &p; }
		return nullptr;
	}
	inline FMultiPoly* GetPolyMut(int Z) noexcept {
		for(auto& p : multi_poly) { if(p.Z == Z) return &p; }
		return nullptr;
	}

	inline std::vector < 
		std::tuple<i32, std::array<double, 2>>
	> GetReferentMeasurements(double x) const {
		std::vector<std::tuple<i32, std::array<double, 2>>> r;
		for(const auto& poly : multi_poly) {
			i32 Z = poly.Z;
			double e_avg_measured = poly::Eval(x, poly.pol);
			double gain = NominalE(Z) / e_avg_measured;
			r.emplace_back(Z, std::array{e_avg_measured, gain});
		}
		return r;
	}

	inline bool IsSaneZ() const {
		return std::is_sorted(multi_poly.begin(), multi_poly.end(), [](const auto& lhs, const auto& rhs) { return lhs.Z < rhs.Z; });	
	}
	inline bool IsSane(double x, std::vector<double>& evals) const {
		/* Check for each 'x' that the each consecutive nominal energy evaluation along the `multi_poly` yields decreasing result.
		 * AKA: that for some `x`, the nominal value of Z cannot be above the value of Z-1. */
		evals.clear();
		for(const auto& mp : multi_poly) {
			evals.push_back(poly::Eval(x, mp.pol));	
		}
		return std::is_sorted(evals.begin(), evals.end());
	}

	enum class Bound { Upper, Lower };
	template<Bound b> 
	auto GetMeasurement(double x, double e) const -> typename Vec::const_iterator {
		static auto comparator = [](const std::pair<double, double> q, const FMultiPoly& p) {
			const auto [x,e] = q;
			return e < poly::Eval(x, p.pol);
		};
		if constexpr(b == Bound::Upper) {
			return std::upper_bound( multi_poly.begin(), multi_poly.end(), std::pair{x,e}, comparator);
		} else {
			return std::lower_bound( multi_poly.begin(), multi_poly.end(), std::pair{x,e}, comparator);
		}
	}

	/* `x` here is the absolute value from [0, 640]. */
	template<ExtrapolateLowZ extrapolate = ExtrapolateLowZ::No>
	double Value(const double x, const double e) const noexcept {
		auto it = this->GetMeasurement<Bound::Upper>(x, e);
		if(it == multi_poly.end()) {
			if(multi_poly.size() == 0) return NAN;
			return FOOTAsicGainParam::NominalE( multi_poly.back().Z ) / poly::Eval(x, multi_poly.back().pol);
		}
		else if(it == multi_poly.begin()) {
			if constexpr(extrapolate == ExtrapolateLowZ::Yes) {
				const auto& [Z1, p1] = *it;
				const auto& [Z2, p2] = *(it+1);
				double e1 = poly::Eval(x, p1);
				double e2 = poly::Eval(x, p2);
				double g1 = NominalE(Z1) / e1;
				double g2 = NominalE(Z2) / e2;
				auto gain_line = GetLine( {e1,g1}, {e2,g2} );
				return poly::Eval(e, gain_line);
			} else {
				return NominalE( multi_poly.front().Z ) / poly::Eval(x, multi_poly.front().pol);
			}
		}
		else {
			const auto& [Z1, p1] = *it;
			const auto& [Z2, p2] = *(it-1);
			double e1 = poly::Eval(x, p1);
			double e2 = poly::Eval(x, p2);
			double g1 = NominalE(Z1) / e1;
			double g2 = NominalE(Z2) / e2;
			auto gain_line = GetLine( {e1,g1}, {e2,g2} );
			return poly::Eval(e, gain_line);
		}
	}
	FOOTAsicGainParam() = default;
	virtual ~FOOTAsicGainParam() = default;
	ClassDef(FOOTAsicGainParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTAsicGainParam, 0);

struct FOOTGainParam {
	static constexpr double PROTON_ADC = FOOTAsicGainParam::PROTON_ADC;
	using AsicArray = std::array<FOOTAsicGainParam, TFOOTMapCont::N_ASIC>;
	GET_HELP_AUX_IMPL;

	ADD_SERIALIZABLE_FIELD(AsicArray, fit, {}, 0);

	/* Throws on invalid access. */ 
	inline const FOOTAsicGainParam& GetASIC(double x) const { return fit.at( x / TFOOTMapCont::N_STRIPS_PER_ASIC ); }

	template<ExtrapolateLowZ extrapolate = ExtrapolateLowZ::No>
	double CorrectionFactor(double x, double e) const noexcept {
		const auto& asic = GetASIC(x);
		return asic.Value<extrapolate>(x, e);
	}

	inline bool IsSane() const {
		std::vector<double> evals(20);
		for(int i=0; i<TFOOTMapCont::N_ASIC; ++i) {
			const auto& asic = fit[i];			
			if(! asic.IsSaneZ() ) return false;

			for(int n=0; n < TFOOTMapCont::N_STRIPS_PER_ASIC; ++n) {
				double x = TFOOTMapCont::N_STRIPS_PER_ASIC*i + n + 0.5; // centre of the strip. 
				if(! asic.IsSane(x, evals) ) return false;
			}
		}
		return true;
	}

	/* Small helper function to plot what we're actually gain matching upon.
	 * Highly inefficient, but it's w/e. */
	[[ nodiscard ]] inline TH2D* GetHisto (
		int nbins_x = 640,
		int nbins_y = 500,
		int lo_y    = 0,
		int hi_y    = 4000
	) const {
		TH2D* h = new TH2D("_hFOOTGainParam", "FOOT Gain Parameter", 
			nbins_x, 0, TFOOTMapCont::N_STRIPS,
			nbins_y, lo_y, hi_y );
		
		for(int ix = 1; ix <= h->GetNbinsX(); ++ix) {
			double x = h->GetXaxis()->GetBinCenter(ix);
			const auto& asic = GetASIC(x);
			
			for(int iy = 1; iy <= h->GetNbinsY(); ++iy) {
				double y = h->GetYaxis()->GetBinCenter(iy);
				double value = asic.Value(x, y);
				h->SetBinContent(ix, iy, value);
			}
		}
		return h; /* Draw via: `h->Draw("SURF1")`  :-)  */
	}
	[[ nodiscard ]] inline std::pair<TGraph*, TGraph*> GetGraph (
		const double x,
		int nbins_y = 500,
		int lo_y    = 0,
		int hi_y    = 4000
	) const {
		TGraph* g = new TGraph(nbins_y);
		const auto& asic = GetASIC(x);
		for(int i=0; i<nbins_y; ++i) {
			double y = lo_y + (i+0.5) * (hi_y - lo_y) / nbins_y; // centre
			double value = asic.Value(x, y);
			g->SetPoint(i, y, value);
		}
		g->SetLineWidth(4);
		g->SetLineColor(kRed + 2);
		auto m = asic.GetReferentMeasurements(x);
		TGraph* gpts = new TGraph(m.size());
		for(int i=0; i<(int)m.size(); ++i) {
			const auto& [e_ref, gain_ref] = std::get<1>(m[i]);
			gpts->SetPoint(i, e_ref, gain_ref);
		}
		gpts->SetMarkerStyle(20);
		gpts->SetMarkerSize(1.35);
		return { gpts, g };
		/* Draw via:
		 * auto [pts, graph] = GetGraph();
		 * pts->Draw("AP"); graph->Draw("L SAME"); */
	}

	FOOTGainParam() = default;
	virtual ~FOOTGainParam() = default;
	ClassDef(FOOTGainParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(FOOTGainParam, 0);

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
	template<> struct tuple_size<::RNFOOTCluster> : integral_constant<size_t, 5> {};
	template<> struct tuple_element<0, ::RNFOOTCluster> { using type = double; };
	template<> struct tuple_element<1, ::RNFOOTCluster> { using type = double; };
    template<> struct tuple_element<2, ::RNFOOTCluster> { using type = u32; };
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
