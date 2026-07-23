#pragma once

#include "monad/monad.hxx"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>

#include "TFOOTMapCont.h"
#include "util/json_struct_def.hh"
#include "util/PolyFitter.h"
#include "util/Tracking.h"

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
	
	/* Ordered vector, in `Z`. If not, analysis is fucked. Soothe paranoia w/ `IsSane()` method. */
	ADD_SERIALIZABLE_FIELD(Vec, multi_poly, {}, 0);

	inline const FMultiPoly* GetPoly(int Z) const noexcept {
		for(const auto& p : multi_poly) { if(p.Z == Z) return &p; }
		return nullptr;
	}
	inline FMultiPoly* GetPoly(int Z) noexcept {
		for(auto& p : multi_poly) { if(p.Z == Z) return &p; }
		return nullptr;
	}

	inline double GetReferentMeasurement(int Z, double x) const { // slow
		const auto* p = GetPoly(Z);
		if(!p) return NAN;
		return poly::Eval(x, p->pol);
	}

	inline bool IsSaneZ() const {
		return std::is_sorted( multi_poly.begin(), multi_poly.end(), 
			[](const auto& lhs, const auto& rhs) { return lhs.Z < rhs.Z; }
		);	
	}

	inline std::vector<int> GetZ() const {
		std::vector<int> Zs{}; Zs.reserve(multi_poly.size());
		std::transform(multi_poly.begin(), multi_poly.end(), std::back_inserter(Zs), [](const auto& x) { return x.Z; });
		return Zs;
	}

	inline bool IsSane(double x, std::vector<double>& evals) const {
		/* Check for each 'x' that the each consecutive nominal energy evaluation along the `multi_poly` yields decreasing result.
		 * AKA: that for some `x`, the nominal value of Z cannot be above the value of Z-1. */
		evals.clear();
		for(const auto& mp : multi_poly) {
			evals.push_back(poly::Eval(x, mp.pol));	
		}
		return evals.size() > 0 and std::is_sorted(evals.begin(), evals.end());
	}

	/* Return an iterator to the first referent measurment line (along E; for set x) that has ref >= e */
	inline auto GetMeasurementBound(double x, double e) const -> typename Vec::const_iterator {
		return std::lower_bound( multi_poly.begin(), multi_poly.end(), std::pair{x,e}, 
			[](const FMultiPoly& p, const std::pair<double, double>& q) noexcept {
				const auto [x, e] = q;
				return poly::Eval(x, p.pol) < e;
			}
		);
	}
	
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

	inline std::vector<int> GetNominalZ() const {
		std::vector<int> Zs{}; Zs.reserve(nominal_value.size());
		std::transform(nominal_value.begin(), nominal_value.end(), std::back_inserter(Zs), [](const auto& x) { return x.Z; });
		return Zs;
	}

	inline bool IsSane() const {
		std::vector<double> evals(20);
		const auto nominalZ = GetNominalZ();

		for(int i=0; i<TFOOTMapCont::N_ASIC; ++i) {
			const auto& asic = fit[i];			
			if(!asic.IsSaneZ()) return false;
			
			const auto Zs = asic.GetZ();
			if(Zs != nominalZ) return false;

			for(int n=0; n < TFOOTMapCont::N_STRIPS_PER_ASIC; ++n) {
				double x = TFOOTMapCont::N_STRIPS_PER_ASIC*i + n + 0.5; // centre of the strip. 
				if(!asic.IsSane(x, evals)) return false;
			}
		}
		/* Check also if the referent Z measurements are sorted, and all are unique. */	
		return ( 
			std::is_sorted(nominal_value.begin(), nominal_value.end()) and
			std::adjacent_find(nominal_value.begin(), nominal_value.end(), 
				[](const auto& lo, const auto& hi) { return lo.Z == hi.Z; }) == nominal_value.end()
		);
	}

	inline double GetNominalValue(i32 Z) const noexcept {
		auto it = std::find_if( nominal_value.begin(), nominal_value.end(), [Z](const auto& x) { return x.Z == Z; });	
		if(it == nominal_value.end()) return NAN;
		return it->val;
	}

	double CorrectionFactor(double x, double e) const noexcept {
		const auto& asic = GetASIC(x);
		double gain;
		
		auto ref_hi = asic.GetMeasurementBound(x, e);
		if(ref_hi == asic.multi_poly.begin()) { // `e` smaller than the lower bound.
			const auto& [Z, poly] = *ref_hi;
			double avg_val = GetNominalValue(Z);
			double val_x = poly::Eval(x, poly);
			gain = avg_val / val_x;
		} else if(ref_hi == asic.multi_poly.end()) { // `e` bigger than the upper bound.
			const auto& [Z, poly] = asic.multi_poly.back();
			double avg_val = GetNominalValue(Z);
			double val_x = poly::Eval(x, poly);
			gain = avg_val / val_x;
		}
		else { // Iterator valid and not pointed at the range head.
			const auto& [Z_lo, poly_lo] = *(ref_hi - 1);
			const auto& [Z_hi, poly_hi] = *ref_hi;
			
			double nom_lo = GetNominalValue(Z_lo);
			double nom_hi = GetNominalValue(Z_hi);
			double evaluated_lo = poly::Eval(x, poly_lo);
			double evaluated_hi = poly::Eval(x, poly_hi);
			double gain_lo = nom_lo / evaluated_lo;
			double gain_hi = nom_hi / evaluated_hi;
			auto line = GetLine2D( {evaluated_lo, gain_lo}, {evaluated_hi, gain_hi} );
			gain = poly::Eval(e, line.array());
		}
		return gain;
	}

	/* Small helper function to plot what we're actually gain matching upon.
	 * Highly inefficient, but it's w/e. */
	[[ nodiscard ]] inline TH2D* GetHisto (
		int nbins_x = 640,
		int nbins_y = 500,
		int lo_y    = 0,
		int hi_y    = 2500
	) const {
		TH2D* h = new TH2D("_hFOOTGainParam", "FOOT Gain Parameter", 
			nbins_x, 0, TFOOTMapCont::N_STRIPS,
			nbins_y, lo_y, hi_y );
		
		for(int ix = 1; ix <= h->GetNbinsX(); ++ix) {
			double x = h->GetXaxis()->GetBinCenter(ix);
			
			for(int iy = 1; iy <= h->GetNbinsY(); ++iy) {
				double y = h->GetYaxis()->GetBinCenter(iy);

				double value = this->CorrectionFactor(x, y);
				h->SetBinContent(ix, iy, value);
			}
		}
		h->GetXaxis()->SetTitle("Strip #");
		h->GetYaxis()->SetTitle("Cluster ADC");

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
			double value = this->CorrectionFactor(x, y);
			g->SetPoint(i, y, value);
		}
		g->SetLineWidth(4);
		g->SetLineColor(kRed + 2);
		const std::vector<i32> Zs = this->GetNominalZ();
		TGraph* gpts = new TGraph{};
		for(i32 Z : Zs) {
			double value = asic.GetReferentMeasurement(Z, x);
			double gain = this->CorrectionFactor(x, value);
			gpts->AddPoint(value, gain);
		}

		gpts->SetMarkerStyle(20);
		gpts->SetMarkerSize(1.35);
		return { gpts, g };
		/* Draw via:
		 * auto [pts, graph] = GetGraph();
		 * graph->Draw("AL"); pts->Draw("P SAME"), ; */
	}
	[[ nodiscard ]] inline std::pair<TGraph*, TLine*> GetRefZGraph(
		int Z, 
		int Npts = 640
	) const {
		TGraph* g = new TGraph(Npts);
		double nm = this->GetNominalValue(Z);
		if(!std::isfinite(nm)) 
			throw std::invalid_argument(Form("FOOTGainParam::GetRefZGraph, Passed Z=%d, which isn't found in the config.", Z));
		double xlo = 0, xhi = 640;
		TLine *lref = new TLine(xlo, nm, xhi, nm);
		for(int i=0; i<Npts; ++i) {
			double x = xlo + (i+0.5)*(xhi - xlo) / Npts;
			const auto& asic = GetASIC(x);
			double val = asic.GetReferentMeasurement(Z, x);
			g->SetPoint(i, x, val);
		}
		constexpr double lwidth_ = 4.0;
		g->SetLineWidth(lwidth_);
		g->SetLineColor(kRed + 2);
		lref->SetLineStyle(2);
		lref->SetLineColor(kMagenta + 2);
		lref->SetLineWidth(lwidth_);

		return {g, lref};
	}
	[[ nodiscard ]] inline std::vector<std::pair<TGraph*, TLine*>> GetAllRefZGraph(
		int Npts = 640
	) const {
		std::vector<std::pair<TGraph*, TLine*>> r;
		const auto nominalZ = GetNominalZ();
		for(const i32 Z : nominalZ) 
			r.emplace_back( GetRefZGraph(Z, Npts) );
		return r;
	}

	FOOTGainParam() = default;
	virtual ~FOOTGainParam() = default;
	ClassDef(FOOTGainParam, 1);
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

	inline Orientation GetOrientation() const noexcept {
		if(orientation[0] == 'x' or (orientation.length() > 1 and orientation[1] == 'x')) return Orientation::X;
		if(orientation[0] == 'y' or (orientation.length() > 1 and orientation[1] == 'y')) return Orientation::Y;
		return Orientation::UNKNOWN;
	};

	inline double R() const noexcept {
		if(orientation[0] == '-') return -1.0;
		return 1.0;
	}

	/* Convert raw cluster position `cx` and its integrated ADC value
	 * `ce` into the final calibrated ADC value. */
	double E(const RNFOOTCluster& ) const noexcept;

	/* Main method: convert cluster energy (E) to nominal charge (Q)
	 * If the dependence is E(Q) = A * Q^a, then:
	 * f = 1/A, c = 1/a <=> Q(E) = (f * E)^c 
	 * This energy, has to be properly both gain-matched and delta-corrected. */
	inline double Q(double E) const noexcept {
		if(!this->is_initialized_)
			QParamInit();
		return std::pow(f_*E, c_);
	}
	double Q(const RNFOOTCluster& ) const noexcept;

	/* Calculate the (bare) cluster hit position in millimeters. Note: this does not 
	 * take offsets or rotations into account. */
	double BarePosition(const RNFOOTCluster& ) const noexcept;

	/* Calculate the cluster hit position in millimeters, along the detector's axis. */
	double X0(const RNFOOTCluster& ) const noexcept;

	/* Get params for the nominal dependence: E(Q) = <0> * Q ^ <1> */
	template<size_t N>
	double GetQParam() const noexcept {
		if(!this->is_initialized_)
			QParamInit();
		if constexpr(N == 0) {
			return 1.0 / f_; 
		} else if constexpr(N == 1) {
			return 1.0 / c_;
		} else {
			static_assert(N < 2, "Template param must be 0 or 1.");
		}
	}
	inline void ResetQ() const noexcept { this->is_initialized_ = false; }
	
protected: 
	/* Some cached values for quick Q- calculation.
	 * NB: if the object is re-evaluted, the values *need* to be recomputed, but 
	 * JSON propagator cannot know this. Meaning that `ResetQ` has to be called manually. */
	mutable double f_, c_;
	mutable bool is_initialized_ = 0;

	inline void QParamInit() const {
		std::vector<double> x, y;
		for(auto [Q,E] : gain.nominal_value) {
			x.push_back( std::log(Q) );
			y.push_back( std::log(E) );
		}
		auto r = PolyFit<1>(x,y);
		this->f_ = std::exp(-r[0]);
		this->c_ = 1.0 / r[1];
		this->is_initialized_ = true;
	}

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
	bool IsOk() const noexcept { return std::isfinite(delta); }

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
	};

	double fCX = 0; /* Cluster mean strip position. */
	double fCE = 0; /* Cluster summed energy. */
	u32    fCM = 0; /* Cluster multiplicity. */
	ClusterType fCT{}; /* Cluster type. */
	
	FOOTClusterFit fit{};

	inline double Delta() const noexcept { return mnd::rround<int>(fCX) - fCX; }
 
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
    template<> struct tuple_element<2, ::RNFOOTCluster> { using type = u32; };
    template<> struct tuple_element<3, ::RNFOOTCluster> { using type = RNFOOTCluster::ClusterType; };
}

struct RNFOOTCal {
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

inline void Add(FOOTParam&, const FOOTParam&) {}
inline void Add(FOOTBoxParam&, const FOOTBoxParam&) {}

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
