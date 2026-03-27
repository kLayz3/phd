#pragma once

#include "TAxis.h"
#include "TH2D.h"
#include <array>

struct Point { double x, y; };

inline std::array<double, 2> GetLine(const Point& p1, const Point& p2) noexcept {
	double slope = (p2.y - p1.y) / (p2.x - p1.x);
	double offset = -slope * p1.x + p1.y;
	return { offset, slope };
}
inline std::array<double, 2> GetLine(const std::array<std::array<double, 2>, 2> p12) noexcept {
	return GetLine( {.x = p12[0][0], .y = p12[0][1]}, {.x = p12[1][0], .y = p12[1][1]} );
}

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

