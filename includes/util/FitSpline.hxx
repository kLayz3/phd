#pragma once

#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "PolyFitter.h"

#include "json_struct_def.hh"
#include <iostream>


template<std::size_t R>
std::array<double, R+1> FitSpline(TH2D* h2) {
	auto pfx = std::unique_ptr<TProfile>( h2->ProfileX("pfx__") );
	if(!pfx) {
		throw std::runtime_error("FitSpline(): ProfileX() failed");
	}
	pfx->SetDirectory(nullptr);

	const int N = pfx->GetNbinsX();
	
	std::vector<double> x(N);
	std::vector<double> y(N);
	for(int i=1; i <= N; ++i) {
		x.at(i-1) = pfx->GetBinCenter(i) ;
		y.at(i-1) = pfx->GetBinContent(i) ;
	}
	
	return PolyFit<R>(x,y);
}

template<std::size_t R>
[[ nodiscard ]]
std::pair<
	std::array<double, R+1>,
	TGraph*
> FitSplineAndGraph(TH2D* h2, int Npts = 80) {
	assert(Npts > 0);
	auto r = FitSpline<R>(h2);

	const double xmin = h2->GetXaxis()->GetXmin();
	const double xmax = h2->GetXaxis()->GetXmax();
	
	TGraph* g = new TGraph(Npts);
	for(int i=0; i<Npts; ++i) {
		double x0 = xmin + (i+0.5) * (xmax - xmin) / Npts; // centre
		double y0 = poly::Eval(x0, r);
		g->SetPoint(i, x0, y0);
	}
	
	g->SetLineColor(kRed);
	g->SetLineWidth(3);

	return { r, g };
}
