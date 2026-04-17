#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/FitSpline.hxx"
#include "../../includes/util/Tracking.hxx"
#include "../../includes/util/MacroHelpers.hxx"

constexpr const char* label[] = {
	"21", "22", "31", "41"
}; constexpr int NSCI = static_cast<int>( sizeof(label)/sizeof(*label) );
enum class DoSave { no, yes };

void sci_x_cal (
	std::string fileName = "", 
	int i_sci = 0,
	A3 binning_x = {100, -30, 30},
	A3 binning_d = {100, -10, 10},
	A2 fit_range = {-5, 5},
	A2 sci21_cut = {-DBL_MAX, DBL_MAX},
	A2 sci22_cut = {-DBL_MAX, DBL_MAX},
	A2 sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no
) {
	if(i_sci != 0 and i_sci != 1)
		throw std::invalid_argument("Second arg (i_sci) must be either 0 or 1. Calibrating only SCI21 and SCI22 here.");

	ROOT::EnableImplicitMT();
	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_params = f->Get < 
			std::remove_reference_t<decltype(*tpc_params)>
		> ("FRS_tpc_parameters");
		if(!tpc_params)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
		sci_params = f->Get < 
			std::remove_reference_t<decltype(*sci_params)>
		> ("FRS_sci_parameters");
		if(!sci_params)
			throw std::runtime_error(Form("SCI param is nullptr. Fix it (line: %d).", __LINE__));

	}
	const auto& sci_param = sci_params->at(i_sci); 
	const double z0 = sci_param.z0; 
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	/* Take TPC21,22,23,24 for the extrapolations. */
	constexpr int N = 4;
	const std::array<double, N> zTPC = {
		tpc_params->at(0).z0,
		tpc_params->at(1).z0,
		tpc_params->at(2).z0,
		tpc_params->at(3).z0,
	};
	/* Optionally skip some TPC's. */
//#define SKIP_TPC21
#define SKIP_TPC22
//#define SKIP_TPC23
#define SKIP_TPC24
	
	std::vector<double> xs, ys, zs;
	xs.reserve(N); ys.reserve(N); zs.reserve(N);
	std::array<double, 2> fx, fy;

	TH2P* hist  = new TH2P(Form("SCI%s X[mm]:TPC extr. [mm]", label[i_sci]), 
		binning_x[0], binning_x[1], binning_x[2],
		binning_x[0], binning_x[1], binning_x[2]); 
	TH2P* histd = new TH2P(Form("SCI%s X -  TPC extr.[mm]:TPC extr. [mm]", label[i_sci]), 
		binning_x[0], binning_x[1], binning_x[2],
		binning_d[0], binning_d[1], binning_d[2]);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);

	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_xy = new TH2P(Form("Referent Y-position [mm]:Referent X-position [mm]@SCI%s", label[i_sci]), 100, -40, 40, 100, -40, 40);
	auto* h2_ab = new TH2P("TPC-derived Y-angle [mrad]:X-angle [mrad]@S2", 100, -20, 20, 100, -20, 20);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		xs.clear(); ys.clear(); zs.clear();
		
		const auto& sci = frs->sci[i_sci];
		if(sci.hits.size() == 0) continue;
		
		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];
		if(sci21.hits.size() != 1) continue;
		if(sci22.hits.size() != 1) continue;
		if(sci31.hits.size() != 1) continue;

		h1_sci21->Fill(sci21.E);
		h1_sci22->Fill(sci22.E);
		h1_sci31->Fill(sci31.E);

		if(!mnd::IsInside(sci21.E, sci21_cut)) continue;
		if(!mnd::IsInside(sci22.E, sci22_cut)) continue;
		if(!mnd::IsInside(sci31.E, sci31_cut)) continue;

		h1_sci21_cut->Fill(sci21.E);
		h1_sci22_cut->Fill(sci22.E);
		h1_sci31_cut->Fill(sci31.E);

		for(int i = 0; i < N; ++i) {
#ifdef SKIP_TPC21
			if(i == 0) continue;
#endif
#ifdef SKIP_TPC22
			if(i == 1) continue;
#endif
#ifdef SKIP_TPC23
			if(i == 2) continue;
#endif
#ifdef SKIP_TPC24
			if(i == 3) continue;
#endif
			const auto& tpc = frs->tpc[i];
			if(tpc.hits[0].size() != 1 or tpc.hits[1].size() != 1) continue;

			double x = tpc.X0();
			double y = tpc.Y0();
			if(!std::isfinite(x) or !std::isfinite(y)) continue;
			
			zs.push_back( zTPC[i] );	
			xs.push_back(x);
			ys.push_back(y);
		}
		if(zs.size() < 2) continue;
		PolyFit<1>(zs, xs, fx);
		PolyFit<1>(zs, ys, fy);
		
		double x_extr = fx[0] + fx[1] * z0;
		double y_extr = fy[0] + fy[1] * z0;
		
		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_xy->Fill(x_extr, y_extr);
		h2_ab->Fill(fx[1]*1000, fy[1]*1000);

		hist ->Fill(x_extr, sci.hits[0].x);
		histd->Fill(x_extr, sci.hits[0].x - x_extr);
	}
	
	auto [rg, gerr, g] = FitSplineAndGraph<1, fit_info::GAUSS_MAX> (
		*hist, fit_range[0], fit_range[1], 40, 3.0 /*, Verbosity::CHATTY */
	);

	auto [offset, slope]  = rg;
	double offset0 = sci_param.x_offset; // already inputted. 
	double slope0 = sci_param.x_factor;   // already inputted. 
	
	char text0[20] = {'\0'};
	sprintf(text0, "For SCI%s:", label[i_sci]);

	char text1[1024] = {'\0'};
	sprintf(text1, "Calculated: graph offset = %.5f, graph slope = %.5f.", offset, slope);
	WARN("%s %s\n", text0, text1);

	char text2[1024] = {'\0'};
	sprintf(text2, "Currently in the setup file: (offset, slope) = (%.5f, %.5f)", sci_param.x_offset, sci_param.x_factor);
	WARN("%s\n", text2);

	double offset1 = (sci_param.x_offset - offset) / slope; // values to be written into file. 
	double slope1 = sci_param.x_factor / slope;             // values to be written into file. 
	char result_slope[1024] = {'\0'};
	char result_offset[1024] = {'\0'};
	sprintf(result_offset, "x_offset: " KCYN "%.5f [mm]" KNRM, offset1);
	sprintf(result_slope, "x_factor: " KCYN "%.5f [mm/tdc]" KNRM, slope1);
	WARN("Recommended: %s\n%s\n", result_offset, result_slope);
	WARN("~~~ This is a relative change of: %.3f%% and %.3f%%\n", 
		100*std::abs((slope1-slope0)/slope0), 100*std::abs((offset1 - offset0)/offset0));

	TCanvas *c = new TCanvas("SCI-pos", "c", 2400, 1400);
	c->Divide(2,2);
	c->cd(1); hist->Draw("COLZ");
	g->Draw("L SAME");
	gerr->Draw("P SAME");

	c->cd(2); histd->Draw("COLZ");
	c->cd(3);
	PLatex(0.06, text0, text1, text2, "Recommended:", result_offset, result_slope);

	std::string s = std::string{"Extrapolation from: TPC"} +
#ifndef SKIP_TPC21
		"21 " +
#endif
#ifndef SKIP_TPC22
		"22 " +
#endif
#ifndef SKIP_TPC23
		"23 " +
#endif
#ifndef SKIP_TPC24
		"24 " +
#endif
		+ "";

	c->cd(4);
	PLatex(0.06, 
		std::move(s), 
#ifndef SKIP_TPC21
		Form("TPC21: @%.1f mm", tpc_params->at(0).z0),
#endif
#ifndef SKIP_TPC22
		Form("TPC22: @%.1f mm", tpc_params->at(1).z0),
#endif
#ifndef SKIP_TPC23
		Form("TPC23: @%.1f mm", tpc_params->at(2).z0),
#endif
#ifndef SKIP_TPC24
		Form("TPC24: @%.1f mm", tpc_params->at(3).z0),
#endif
		Form("SCI%s: @%.1f mm", label[i_sci], sci_param.z0)
	);
	TCanvas* cs = new TCanvas("SCIe", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	
	TCanvas* cTr = new TCanvas("TPC-tracks", "TPC-tracks", 2000, 1200);
	cTr->Divide(2,2);
	cTr->cd(1); h2_track_x->Draw("COLZ");
	cTr->cd(3); h2_track_y->Draw("COLZ");
	cTr->cd(2); h2_xy->Draw("COLZ");
	cTr->cd(4); h2_ab->Draw("COLZ");
	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("SCI%s", label[i_sci]) });
	}
}
