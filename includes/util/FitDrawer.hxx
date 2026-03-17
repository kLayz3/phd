/* PolyFitter is meant to not depend on ROOT,
 * but here is a small wrapper to draw fit results. 
 * Diagnostic small lib. */

#pragma once

#include "PolyFitter.h"
#include <algorithm>
#include <tuple>
#include "TGraph.h"

template<std::size_t R>
[[ nodiscard ]] 
std::pair<TGraph*, TGraph*> FitAndDraw (
	const std::vector<double>& x, 
	const std::vector<double>& y, 
	std::array<double, R+1>& result,
	double ratio_outside = 0.1,
	const int Npts = 60) {
	
	PolyFit<R>(x, y, result);

	TGraph* g0 = new TGraph(x.size(), x.data(), y.data());
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
}

/* First one should be drawn as Draw("P SAME"),
 * second one as a line: Draw("L SAME") */
