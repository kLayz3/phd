#pragma once

#include "TAxis.h"
#include "TH2D.h"
#include <array>

inline void FillTrack(TH2D* hist, const std::array<double, 2>& a) {
	double xmin = hist->GetXaxis()->GetXmin();
	double xmax = hist->GetXaxis()->GetXmax();	

	int NBins = hist->GetNbinsX();
	TAxis* xaxis = hist->GetXaxis();

	double x, y;
	for(int i=1; i <= NBins; ++i) {
		x = xaxis->GetBinCenter(i);
		y = a[0] + a[1]*x;
		hist->Fill(x,y);
	}
}

