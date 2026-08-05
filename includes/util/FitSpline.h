#pragma once

#include "Rtypes.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "PolyFitter.h"
#include "json_struct_def.hh"

#include "GaussFitMax.hxx"
#include <climits>

/* By default, profile will sample the projected bin's mean.
 * But sometimes distributions can be skewed and we just wish to make
 * a quick dirty Gaussian fit around the max value. In this case
 * the second template parameter can be supplied, but this will make
 * the procedure rather hefty as frequent (re)alloc's have to happen. */
namespace mnd::detail {

/* For a TH2D* histogram, get all projections along Y (onto X). A simple
 * case is just a TProfile, which maps each projection to a number (corresponding to the MAX bin).
 *
 * More complex is the variant fit_info::GAUSS_MAX where each projection is treated as its own TH1D,
 * and Gauss is fitted onto that. What gets returned back is a sequence of (x,y,w) where `x` are the bin-centers,
 * the `y` their values of the slice, and `w` weight. */
template < 
	fit_info sliceType = fit_info::PROFILE_MAX
> 
std::tuple <
	std::vector<double>, // x
	std::vector<double>, // y
	std::vector<double>, // w
	TGraphErrors*        // [pts]
> FitSplineImpl (
	TH2D* h2, 
	double& x_lo,
	double& x_hi,
	double side_ratio,
	uint32_t niter,
	Verbosity v
) {
	std::unique_ptr<TH1D> pfx;
	TAxis* xax = h2->GetXaxis();
	x_lo = std::max( xax->GetXmin(), x_lo );
	x_hi = std::min( xax->GetXmax(), x_hi );
	xax->SetRangeUser(x_lo, x_hi);

	if constexpr(sliceType == fit_info::PROFILE_MAX) {
		(void)side_ratio; (void)niter;
		pfx = std::unique_ptr<TProfile>( h2->ProfileX("pfx__") );
		if(!pfx) {
			throw std::runtime_error("FitSpline(): ProfileX() failed");
		}
		pfx->SetDirectory(nullptr);
	}
	else {
		int ibin_first = xax->GetFirst();
		int ibin_last = xax->GetLast();
		double xfirst = xax->GetBinLowEdge(ibin_first);
		double xlast = xax->GetBinUpEdge(ibin_last);
		int nbins = ibin_last - ibin_first + 1;
		
		if(v > 1) { 
			fprintf(stderr, "[%s]: from \'%s\'; taking bin indices: [%d,%d], "
				"with lo = %.1f, hi = %.1f, N=%d\n", 
				"FitSpline", h2->GetTitle(), ibin_first, ibin_last, xfirst, xlast, nbins);
		}

		pfx = std::make_unique<TH1D>("pfx__", "pfx__", 
			nbins, xfirst, xlast);
		pfx->SetDirectory(nullptr);
		for(int ix = ibin_first, i=1; ix <= ibin_last; ++ix, ++i) {
			auto py = std::unique_ptr<TH1D>(h2->ProjectionY(Form("_py__%d", ix), ix, ix));
			py->SetDirectory(nullptr);
			auto [result, err] = GaussFitMax(py.get(), side_ratio, niter, v);
			pfx->SetBinContent(i, result[1]);
			pfx->SetBinError(i, err[1]);
			if(v > 1) {
				fprintf(stderr, "[%d]: %.2f +- %.2f\n", i, result[1], err[1]);
			}
		}
	}
	const int N = pfx->GetNbinsX();
	
	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> w; w.reserve(N);
	for(int i=1; i <= N; ++i) {
		if(std::isfinite( pfx->GetBinContent(i) )) {
			x.push_back( pfx->GetBinCenter(i) );
			y.push_back( pfx->GetBinContent(i));
			double stddev = pfx->GetBinError(i);
			w.push_back( 1.0 / (stddev*stddev) );
		}
	}
	TGraphErrors* gerr = new TGraphErrors( pfx.get() ); 
	gerr->SetMarkerStyle(20);
	gerr->SetMarkerSize(1.2);
	xax->SetRange();
	if(v > 0) {
		std::cerr << "[FitSpline]" << h2->GetName() << '(' 
			<< ((sliceType == fit_info::PROFILE_MAX)? "PROFILE_MAX": "GAUSS_MAX") << ") "
			<< "fitting with " << std::endl
			<< "x: " << x << std::endl
			<< "y: " << y << std::endl
			<< "w: " << w << std::endl;
	}
	return { std::move(x), std::move(y), std::move(w), gerr };
}

/* Once the (x,y,w) sequence was fitted with an Rth degree poly, we want to evaluate this poly
 * and get back a fitted graph. */
template <
	fit_info sliceType, 
	typename Cont
> TGraph* GetGraphImpl(
	const Cont& r, // 
	TH2D* h2,
	double x_lo,
	double x_hi,
	int Npts,
	Verbosity v
) {
	assert((Npts > 0) && "Cannot pass 0 or negative number for points here.");
	if(v > 0) std::cerr << "Result: " << r << std::endl;

	TGraph* g = new TGraph(Npts);
	if(r.size() > 0) {
		for(int i=0; i<Npts; ++i) {
			double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts; // centre
			double y0 = poly::Eval(x0, r);
			if(v > 0) fprintf(stderr, "[%s:%s] Fit: (x0, y0): {%.2f, %.2f}\n", "FitSpline", h2->GetName(), x0, y0);
			g->SetPoint(i, x0, y0);
		}
		if constexpr(sliceType == fit_info::PROFILE_MAX) {
			g->SetLineColor(pCol_);
		} else { 
			g->SetLineColor(gCol_);
		}

		g->SetLineWidth(4);
	}
	return g;
}

} /* namespace mnd::detail */ 

/* Compile time version; poly degree `R` known at compile time. */
template < 
	std::size_t R, 
	fit_info sliceType = fit_info::PROFILE_MAX
> [[ nodiscard ]]
std::tuple <
	std::array<double, R+1>,
	TGraphErrors*, // 'Measured' points
	TGraph*        // Polynomial fit graph
> FitSpline (
	TH2D* h2, 
	double x_lo = -DBL_MAX,
	double x_hi =  DBL_MAX,
	int Npts = 80,
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	uint32_t niter = 2,
	Verbosity v = Verbosity::SILENT
) {
	auto [x,y,w,gerr] = mnd::detail::FitSplineImpl<sliceType>(h2, x_lo, x_hi, side_ratio, niter, v);	
	std::array<double, R+1> r = {NAN};
	if(x.size() >= R+1) PolyFit<R>(x,y,w,r);
	
	TGraph* g = mnd::detail::GetGraphImpl<sliceType>(r, h2, x_lo, x_hi, Npts, v);
	return { r, gerr, g };
}

/* Runtime version, `R` not known at compile time. */
template <
	fit_info sliceType = fit_info::PROFILE_MAX
> [[ nodiscard ]]
std::tuple <
	std::vector<double>,
	TGraphErrors*,
	TGraph*
> FitSpline (
	std::size_t R,
	TH2D* h2, 
	double x_lo = -DBL_MAX,
	double x_hi =  DBL_MAX,
	int Npts = 80,
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	uint32_t niter = 2,
	Verbosity v = Verbosity::SILENT
) {
	auto [x,y,w,gerr] = mnd::detail::FitSplineImpl<sliceType>(h2, x_lo, x_hi, side_ratio, niter, v);	
	std::vector<double> r = {};
	if(x.size() >= R+1) PolyFit(R,x,y,w,r);

	TGraph* g = mnd::detail::GetGraphImpl<sliceType>(r, h2, x_lo, x_hi, Npts, v);
	return { r, gerr, g };
}

extern template std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, TGraphErrors*> 
mnd::detail::FitSplineImpl<fit_info::PROFILE_MAX> (TH2D* , double& , double& , double , uint32_t , Verbosity );

extern template std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, TGraphErrors*> 
mnd::detail::FitSplineImpl<fit_info::GAUSS_MAX> (TH2D* , double& , double& , double , uint32_t , Verbosity );

extern template std::tuple<std::vector<double>, TGraphErrors*, TGraph*> 
FitSpline<fit_info::PROFILE_MAX>(std::size_t, TH2D* , double , double , int , double , uint32_t , Verbosity );

extern template std::tuple<std::vector<double>, TGraphErrors*, TGraph*> 
FitSpline<fit_info::GAUSS_MAX>(std::size_t, TH2D* , double , double , int , double , uint32_t , Verbosity );
