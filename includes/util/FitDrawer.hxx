/* PolyFitter is meant to not depend on ROOT,
 * but here is a small wrapper to return the graphs of
 * actual fit results. Diagnostic small lib. No need to optimise in,
 * since it shouldn't be called anyway in a tight loop. */

#pragma once

#include "PolyFitter.h"
#include <algorithm>
#include "TGraphErrors.h"

[[ nodiscard ]]
inline std::pair <
    TGraphErrors*, // Points (x,y) that went into the fit.
    TGraph*        // Graph of the fitted polynomial
> FitAndDraw (
    size_t R,
	const mnd::span<const double>& x,
	const mnd::span<const double>& y,
    const mnd::span<const double>& w,
	std::vector<double>& result,
	double ratio_outside = 0.1,
	const int Npts = 60
) {
    /* Weights vector can either be left empty, or must match the size of x,y. */
    if(w.empty()) {
	    PolyFit(R, x, y, result);
    } else { // Inside this call asserted w.size() == x.size()
	    PolyFit(R, x, y, w, result);
    }
    
    /* Maybe in later API, `PolyFit` can throw, we construct the graph here 
     * instead, to not have a possible memory leak. */
    TGraphErrors* g0;
    if(w.empty()) {
	    g0 = new TGraphErrors(x.size(), x.data(), y.data());
    } else {
	    g0 = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, w.data());
    }
	g0->SetMarkerStyle(20);
	g0->SetMarkerSize(1.4);

	const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());

	double xlo = *xmin - ratio_outside * (*xmax - *xmin);
	double xhi = *xmax + ratio_outside * (*xmax - *xmin);

	double dx = (xhi - xlo) / (Npts - 1);
	TGraph* g1 = new TGraph(Npts);

	for(int i=0; i<Npts; ++i) {
		double xp = xlo + dx*i;
		double yp = poly::Eval(xp, result);
		g1->SetPoint(i, xp, yp);
	}
	g1->SetLineColor(kRed);
	g1->SetLineWidth(3);
	g1->SetLineStyle(7);
	
	return {g0, g1};
    /* Can customise the lines further, after the call returns. */
}

/* First one should be drawn as Draw("P SAME"),
 * second one as a line: Draw("L SAME") */
