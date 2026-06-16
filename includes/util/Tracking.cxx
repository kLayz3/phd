#include "Tracking.h"
#include "TAxis.h"
#include "TH2D.h"

void FillTrack (
	TH2D* hist, 
	const std::array<double, 2>& a,
	double z_lo,
	double z_hi
) {
	TAxis* xaxis = hist->GetXaxis();

	const int NBins = xaxis->GetNbins();
	const int bin_lo = std::max(xaxis->FindBin(z_lo), 1);
	const int bin_hi = std::min(xaxis->FindBin(z_hi), NBins);

	double x, y;
	for(int i = bin_lo; i <= bin_hi; ++i) {
		x = xaxis->GetBinCenter(i);
		y = a[0] + a[1]*x;
		hist->Fill(x,y);
	}
}

void FillTrack (
	TH2D* hist, 
	const mnd::geom::Line2D& a,
	double z_lo,
	double z_hi
) {
	FillTrack(hist, a.array(), z_lo, z_hi);
}
