/* PolyFitter is meant to not depend on ROOT,
 * but here is a small wrapper to return the graphs of
 * actual fit results. Diagnostic small lib. No need to optimise in,
 * since it shouldn't be called anyway in a tight loop. */

#pragma once

#include "PolyFitter.h"
#include <algorithm>
#include <string>
#include <stdexcept>
#include "TGraphErrors.h"

using namespace std::literals;

[[ nodiscard ]]
inline std::pair <
    TGraphErrors*, // Points (x,y) that went into the fit.
    TGraph*        // Graph of the fitted polynomial
> FitAndDraw (
    const std::size_t R,
	const mnd::span<const double>& x,
	const mnd::span<const double>& y,
    const mnd::span<const double>& w,
	const mnd::span<const uint32_t>& ignored_indices,
	std::vector<double>& result,
	const double ratio_outside = 0.1,
	const int Npts = 60
) {
	/* Weights vector can either be left empty, or must match the size of x,y. */
	if(!w.empty() && w.size() != x.size())
		throw std::invalid_argument(
			"FitAndDraw(..) => Passed in weights sized"s + std::to_string(w.size()) +
			" which doesn't match points sized " + std::to_string(x.size()));
	if(x.size() != y.size())
		throw std::invalid_argument(
			"FitAndDraw(..) => Passed in x- sized"s + std::to_string(x.size()) +
			" which doesn't match y- sequence sized " + std::to_string(y.size()));

	const uint32_t N = static_cast<uint32_t>(x.size());

	std::unique_ptr<TGraphErrors> g0 {};
    if(w.empty()) {
	    g0 = std::make_unique<TGraphErrors>( x.size(), x.data(), y.data() );
    } else {
	    g0 = std::make_unique<TGraphErrors>( x.size(), x.data(), y.data(), nullptr, w.data() );
    }
	g0->SetMarkerStyle(20);
	g0->SetMarkerSize(1.4);

	const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());

	std::vector<double> x1, y1, w1;
	x1.reserve(N), y1.reserve(N), w1.reserve(N);

	for(uint32_t i=0; i<N; ++i) {
		if( std::find(ignored_indices.begin(), ignored_indices.end(), i)
			== ignored_indices.end()
		) {
			x1.push_back(x[i]);
			y1.push_back(y[i]);
			if(!w.empty())
				w1.push_back(w[i]);
		}
	}
	if(x1.size() < R)
		throw std::invalid_argument(
			"FitAndDraw(..) => After filtering, number of points remaining is"s + std::to_string(x1.size()) +
			" which is less than the polynomial degree to be fitted " + std::to_string(R));

    if(w.empty()) {
	    PolyFit(R, x1, y1, result);
    } else {
	    PolyFit(R, x1, y1, w1, result);
    }

	double xlo = *xmin - ratio_outside * (*xmax - *xmin);
	double xhi = *xmax + ratio_outside * (*xmax - *xmin);

	double dx = (xhi - xlo) / (Npts - 1);
	auto g1 = std::make_unique<TGraph>(Npts);

	for(int i=0; i<Npts; ++i) {
		double xp = xlo + dx*i;
		double yp = poly::Eval(xp, result);
		g1->SetPoint(i, xp, yp);
	}
	g1->SetLineColor(kRed);
	g1->SetLineWidth(3);
	g1->SetLineStyle(7);
	
	return { g0.release(), g1.release() };
}
/* ... can customise the line objects further. */

/* First one should be drawn as Draw("P SAME"),
 * second one as a line: Draw("L SAME") */
