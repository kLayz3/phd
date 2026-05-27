#include "Tracking.h"
#include "TAxis.h"
#include "TH2D.h"

void FillTrack(TH2D* hist, const std::array<double, 2>& a) {
	int NBins = hist->GetNbinsX();
	TAxis* xaxis = hist->GetXaxis();

	double x, y;
	for(int i=1; i <= NBins; ++i) {
		x = xaxis->GetBinCenter(i);
		y = a[0] + a[1]*x;
		hist->Fill(x,y);
	}
}
void FillTrack(TH2D* hist, const mnd::geom::Line2D& a) {
	int NBins = hist->GetNbinsX();
	TAxis* xaxis = hist->GetXaxis();

	double x, y;
	for(int i=1; i <= NBins; ++i) {
		x = xaxis->GetBinCenter(i);
		y = a[0] + a[1]*x;
		hist->Fill(x,y);
	}
}
