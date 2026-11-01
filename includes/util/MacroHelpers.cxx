#include "util/MacroHelpers.h"

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sys/wait.h>
#include <unistd.h>

#include "TH2D.h"
#include "util/Geometry.h"

using PipeDeleter = int (*)(FILE*);

template void canvas::save_all<canvas::Macro>(canvas::Extension , std::vector<std::string_view> );
template void canvas::save_all<canvas::Exe  >(canvas::Extension , std::vector<std::string_view> );
template void canvas::save_all<canvas::Macro>(std::vector<canvas::Extension> , std::vector<std::string_view> );
template void canvas::save_all<canvas::Exe  >(std::vector<canvas::Extension> , std::vector<std::string_view> );

void canvas::DumpPrimitives(TVirtualPad* pad, int depth) {
	if(!pad) return;

	TIter next(pad->GetListOfPrimitives());
	while(TObject* obj = next()) {
		std::cout
			<< std::string(depth * 2, ' ')
			<< obj->ClassName()
			<< "  "
			<< obj->GetName()
			<< "  "
			<< obj->GetTitle()
			<< '\n';

		if(auto* subpad = dynamic_cast<TVirtualPad*>(obj))
			DumpPrimitives(subpad, depth + 1);
	}
}

std::vector<TH1*> canvas::GetHistograms(TVirtualPad* pad) {
	std::vector<TH1*> out {};
	if(!pad) return out;

	TIter next(pad->GetListOfPrimitives());
	while(auto* obj = next()) {
		if(auto* h = dynamic_cast<TH1*>(obj))
			out.push_back(h);

		else if(auto* subpad = dynamic_cast<TVirtualPad*>(obj)) {
			auto sub = GetHistograms(subpad);
			out.insert(out.end(), sub.begin(), sub.end());
		}
	}

	return out;
}

void mnd::FillTrack (
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

void mnd::FillTrack (
	TH2D* hist,
	const mnd::geom::Line2D& a,
	double z_lo,
	double z_hi
) {
	::mnd::FillTrack(hist, a.array(), z_lo, z_hi);
}
