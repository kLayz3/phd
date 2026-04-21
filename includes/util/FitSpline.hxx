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

static const char* module_name_ = "FitSpline";

/* By default, profile will sample the projected bin's mean.
 * But sometimes distributions can be skewed and we just wish to make
 * a quick dirty Gaussian fit around the max value. In this case
 * the second template parameter can be supplied, but this will make
 * the procedure rather hefty as frequent (re)alloc's have to happen. */
template < 
	std::size_t R, 
	fit_info sliceType = fit_info::PROFILE_MAX
> 
std::pair <
	std::array<double, R+1>,
	TGraphErrors*
> FitSpline (
	TH2D* h2, 
	double x_lo = -DBL_MAX,
	double x_hi =  DBL_MAX,
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	Verbosity v = Verbosity::SILENT
) {
	std::unique_ptr<TH1D> pfx;
	TAxis* xax = h2->GetXaxis();
	x_lo = std::max( xax->GetXmin(), x_lo );
	x_hi = std::min( xax->GetXmax(), x_hi );
	xax->SetRangeUser(x_lo, x_hi);

	if constexpr(sliceType == fit_info::PROFILE_MAX) {
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
			printf("[FitSpline] %s: from \'%s\'; taking bin indices: [%d,%d], "
				"with lo = %.1f, hi = %.1f, N=%d\n", 
				module_name_, h2->GetTitle(), ibin_first, ibin_last, xfirst, xlast, nbins);
		}

		pfx = std::make_unique<TH1D>("pfx__", "pfx__", 
			nbins, xfirst, xlast);
		pfx->SetDirectory(nullptr);
		for(int ix = ibin_first, i=1; ix <= ibin_last; ++ix, ++i) {
			auto py = std::unique_ptr<TH1D>(h2->ProjectionY(Form("_py__%d", ix), ix, ix));
			py->SetDirectory(nullptr);
			auto [result, err] = GaussFitMax(py.get(), side_ratio, v);
			pfx->SetBinContent(i, result[1]);
			pfx->SetBinError(i, err[1]);
			if(v > 1) {
				printf("[%d]: %.2f +- %.2f\n", i, result[1], err[1]);
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
		std::cout << "Fitting" << std::endl
			<< "x: " << x << std::endl
			<< "y: " << y << std::endl
			<< "w: " << w << std::endl;
	}
	
	if(x.size() >= R+1) return { PolyFit<R>(x,y,w), gerr };
	//if(x.size() >= R+1) return { PolyFit<R>(x,y), gerr };
	else return { {}, gerr };
}

template <
	std::size_t R,
	fit_info sliceType = fit_info::PROFILE_MAX
> [[ nodiscard ]] 
std::tuple <
	std::array<double, R+1>,
	TGraphErrors*,
	TGraph*
> FitSplineAndGraph (
	TH2D* h2, 
	double x_lo = -DBL_MAX,
	double x_hi =  DBL_MAX,
	int Npts = 80,
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	Verbosity v = Verbosity::SILENT
) {
	assert((Npts > 0) && "Cannot pass 0 or negative number for points here.");
	auto [r, gerr] = FitSpline<R, sliceType>(h2, x_lo, x_hi, side_ratio, v);

	if(v > 0) std::cout << "Result: " << r << std::endl;
	TAxis* xax = h2->GetXaxis();
	x_lo = std::max( xax->GetXmin(), x_lo );
	x_hi = std::min( xax->GetXmax(), x_hi );
	
	TGraph* g = new TGraph(Npts);
	for(int i=0; i<Npts; ++i) {
		double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts; // centre
		double y0 = poly::Eval(x0, r);
		if(v > 0) printf("[%s] Fit: (x0, y0): {%.2f, %.2f}\n", h2->GetName(), x0, y0);
		g->SetPoint(i, x0, y0);
	}
	if constexpr(sliceType == fit_info::PROFILE_MAX) {
		g->SetLineColor(pCol_);
	} else { 
		g->SetLineColor(gCol_);
	}

	g->SetLineWidth(4);

	return { r, gerr, g };
}
