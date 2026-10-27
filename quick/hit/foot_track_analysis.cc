#include "TFRSCalCont.h"
#include "util/CLI.h"
#include "util/Geometry.h"
#include "util/PrettyHisto.h"
#include "util/Tracking.h"
#include "util/MacroHelpers.h"
#include "util/TextUtils.h"

#include "TStyle.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

int main(int argc, char* argv[]) {
	TClass* _cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!_cl or !_cl->GetDataMember("_x"))
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	CLI::App app{"Gate on straight primary beam tracks, and try to estimate the corresponding Cq, Cp, Ct coefficients, \
	              Based on already preliminary results from the Kalman filter. Basically, here we can gate explicitly \
                  on specific \'Z\' tracks (can be wrong z tho)"};
	
	std::string fileName {};
	std::map<u32, A2> cut_q = {};
	A3 binning_q = {120, 0.0, 7.5};
	A3 binning_dx = {100,-5,5};
	A3 binning_dy = {100,-5,5};
	A3 binning_dq = {120, -2.1, 2.1};
	A2 sci21_cut = {NAN,NAN};
	A2 sci22_cut = {NAN,NAN};
	A2 sci31_cut = {NAN,NAN};
	std::vector<canvas::Extension> save = {};
	
	add_logged_option(app, "-f,--file", fileName, "Pass a file name.")
		->check(CLI::ReadPermissions);

	add_logged_option(app, "-q,--cut-q", cut_q, Form("Charge (Z) cut on each FOOT index, index=0,1,..%u.", N_FOOT_DETECTORS));
	add_logged_option(app, "--binning-q", binning_q, "Binning of the charge measurement, Q of each FOOT, in [Z].")
		->delimiter(',');
	add_logged_option(app, "--binning-dx", binning_dx, "Binning of the difference δX (foot - track), in [mm].")
		->delimiter(',');
	add_logged_option(app, "--binning-dy", binning_dy, "Binning of the difference δY (foot - track), in [mm].")
		->delimiter(',');
	add_logged_option(app, "--binning-dq", binning_dq, "Binning of the difference δQ (foot - track), in [Z].")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21", sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22", sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31", sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

#ifdef MND_FOOTTRACK_DEBUG
	if(fileName.empty()) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}
	constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	
	TApplication rootApp("app", 0, 0);
	
	double Cr, Cq, Ct, max_cost, max_cost_f;
	std::array<TH1I*, N_PAIRS> h1_diff_q;
	std::array<TH1I*, N_PAIRS> h1_diff_r;
	std::array<TH1I*, N_PAIRS> h1_diff_t;
	std::array<TH1I*, N_PAIRS> h1_acc_q;
	std::array<TH1I*, N_PAIRS> h1_acc_r;
	std::array<TH1I*, N_PAIRS> h1_acc_t;
	{
		A3* c;
		TParameter<double>* m;
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2);
		max_cost = m->GetVal();
		get_obj(f, m, "FOOT_max_cost_f");
		max_cost_f = m->GetVal();
		for(u32 n=0; n<N_PAIRS; ++n) {
            get_obj(f, h1_diff_r.at(n), Form("FOOT_diff_r_%u", n));
			get_obj(f, h1_diff_q.at(n), Form("FOOT_diff_q_%u", n));
            get_obj(f, h1_diff_t.at(n), Form("FOOT_diff_t_%u", n));
            get_obj(f, h1_acc_r.at( n), Form("FOOT_acc_r_%u", n));
            get_obj(f, h1_acc_q.at( n), Form("FOOT_acc_q_%u", n));
            get_obj(f, h1_acc_t.at( n), Form("FOOT_acc_t_%u", n));
		}
	}
	
	TH1P *resx[N_PAIRS], *resy[N_PAIRS], *resq[N_PAIRS], *h1q[N_PAIRS], *h1_footx[N_PAIRS], *h1_footy[N_PAIRS];
	for(u32 i=0; i<N_PAIRS; ++i) {
		resx[i] = new TH1P(Form("((h1_fitrx_%d))Fit residue [mm]@FOOT%d X", i, i),
			kMagenta-9, binning_dx[0], binning_dx[1], binning_dx[2]);
		resy[i] = new TH1P(Form("((h1_fitry_%d))Fit residue [mm]@FOOT%d Y", i, i),
			kYellow-9, binning_dy[0], binning_dy[1], binning_dy[2]);
		resq[i] = new TH1P(Form("((h1_fitrq_%d))Fit residue [charge unit]@FOOT%d Q", i, i),
			kCyan-9, binning_dq[0], binning_dq[1], binning_dq[2]);

		h1q[i] = new TH1P(Form("((h1q_%d))Q value@FOOT%d", i, i),
			kBlue-9, binning_q[0], binning_q[1], binning_q[2]);
		h1_footx[i] = new TH1P(Form("((h1_footx_%d))FOOT%d X [mm]", i, i),
			kMagenta-9, 320, -50, 50);
		h1_footy[i] = new TH1P(Form("((h1_footy_%d))FOOT%d Y [mm]", i, i),
			kMagenta-9, 320, -50, 50);
	}
	TH2P* h2_resx = new TH2P("((h2_fitrx))Fit residue:FOOT ID@X-orientation",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_dx[0], binning_dx[1], binning_dx[2]);
	TH2P* h2_resy = new TH2P("((h2_fitry))Fit residue:FOOT ID@Y-orientation",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_dy[0], binning_dy[1], binning_dy[2]);
	TH2P* h2_resq = new TH2P("((h2_fitrq))Fit residue:FOOT ID@Charge (Q)",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_dq[0], binning_dq[1], binning_dq[2]);

	TH2P* h2q = new TH2P("((h2q))Q value:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_q[0], binning_q[1], binning_q[2]);

	TH2P* h2_kq = new TH2P("Kq cost [Q^2]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 500, 0, 1);
	TH2P* h2_kr = new TH2P("Kr cost [mm^2]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 400, 0, 10);
	TH2P* h2_kq_sqrt = new TH2P("sqrt_Kq cost [Q]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 500, 0, 1);
	TH2P* h2_kr_sqrt = new TH2P("sqrt_Kr cost [mm]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 400, 0, 3.1628);

	TH2P* h2_ktxy = new TH2P("FOOT - Upstream Y [mm]:FOOT - Upstream x [mm]@target",
		400,-10,10, 400,-10,10);
	TH1P* h1_kt = new TH1P("Kt cost [mm^2]", kYellow-3,
		400, 0, 10);
	TH1P* h1_kt_sqrt = new TH1P("sqrt_Kt cost [mm]", kYellow-3,
		400, 0, 3.162);
	TH1P* h1_diff_upstr_down = new TH1P("Distance Upstream to Downstream track [mm]", kGreen-1,
		400, 0, 100);
	TH1P* h1_has_upstream = new TH1P("Upstream track present [0=false, 1=true]", kYellow-3,
		4, -0.75, 1.25);
	TH1P* h1_score = new TH1P("Track score [a.u.]", kGreen-1,
		500, 0, 50);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", 0xCC70DD_c, 500, 300, 4000);
	auto* h1_sci21_cut  = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut  = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut  = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0xCC80FD_c, 500, 300, 4000);
	auto* h1_sci21_cut2 = new TH1P("((h1_cut2)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut2 = new TH1P("((h1_cut2)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut2 = new TH1P("((h1_cut2)) SCI31 QDC mean [QDC units]@With cut", 0xCC80FD_c, 500, 300, 4000);

	for(const auto& [foot_i, _] : cut_q) {
		if(foot_i >= N_PAIRS)
			ERROR("Supplied index: %u >= %u as Q-cut pair index.\n", foot_i, N_PAIRS);
	}

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTHit>("FOOT");
	auto frs = model->MakeField<RNFRSHit>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h104", fileName);
    const size_t nentries = ntuple->GetNEntries();

	/* Fitting containers, for the track made by the 3 FOOTs projecting onto the 4th. */
	std::array<double, 3> xf, yf, zf, qf;
	ProgressBar bar {
		option::BarWidth{50},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("FOOT track analysis (per event: %s)", fileName.c_str())},
			option::ForegroundColor{Color::red},
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector{FontStyle::bold}}
	};

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		mnd::PrintProgress(bar, entryId, nentries, 500);

		const auto& sci21 = frs->cal.sci[RNFRSCal::SCI21_I];
		const auto& sci22 = frs->cal.sci[RNFRSCal::SCI22_I];
		const auto& sci31 = frs->cal.sci[RNFRSCal::SCI31_I];
		
		if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
		if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
		if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);
		
		/* Skip the event entirely in case a SCI cut isn't met. */
		if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
		if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
		if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

		if(sci21.hits.size() >= 1) h1_sci21_cut->Fill(sci21.E);
		if(sci22.hits.size() >= 1) h1_sci22_cut->Fill(sci22.E);
        if(sci31.hits.size() >= 1) h1_sci31_cut->Fill(sci31.E);

		const double x0_upst = frs->xT;
		const double y0_upst = frs->yT;
		
		const mnd::geom::Line3D lu = RNTrackToLine3D( frs->s2_bt );
		h1_has_upstream->Fill( static_cast<int>(lu.HasValue()) );

		const RNFOOTTrack& t = foot->heavy_fragment;
		if(t.n != N_PAIRS) continue;

		/* Check if charge veto passed. */
		if(!std::all_of(cut_q.begin(), cut_q.end(),
			[&t](const auto& e) -> bool {
				const auto& [foot_i, cut] = e;
				const double qv = t._q[ foot_i ];
				return mnd::IsInside(qv, cut);
			})
		) {
			continue;
		}

		for(size_t i=0; i<N_PAIRS; ++i) {
			double z = t._z[i];
			
			int w=0;
			
			/* Populate the arrays for N-1 fit. */
			for(size_t j=0; j<N_PAIRS; ++j) {
				if(j == i) continue;
				xf[w] = t._x[j];
				yf[w] = t._y[j];
				zf[w] = t._z[j];
				qf[w] = t._q[j];
				++w;
			} if(w != 3) ERROR("w must be 3, hello?");

			auto fx = PolyFit<1>(zf, xf);
			auto fy = PolyFit<1>(zf, yf);
			//auto fx = std::array<double, 2> { t.x0, t.ax };	
			//auto fy = std::array<double, 2> { t.y0, t.ay };

			double x_extr = fx[0] + fx[1]*z;
			double y_extr = fy[0] + fy[1]*z;
			auto [q_extr, track_qvar] = mnd::mean_var(qf);
			
			double dx = t._x[i] - x_extr;
			double dy = t._y[i] - y_extr;
			double dq = t._q[i] - q_extr;
			double sq = t._sq[i]; // is sigma == sqrt(var);

			h1_footx[i]->Fill( t._x[i] );
			h1_footy[i]->Fill( t._y[i] );
			resx[i]->Fill(dx);
			resy[i]->Fill(dy);
			resq[i]->Fill(dq);
			h2_resx->Fill(i, dx);
			h2_resy->Fill(i, dy);
			h2_resq->Fill(i, dq);

			h1q[i]->Fill(t._q[i]);
			h2q->Fill(i, t._q[i]);
			
			/* Find the 3-measurement 'test' track charge params.. */
			h2_kq->Fill( i, dq*dq + sq*sq );
			h2_kr->Fill( i, dx*dx + dy*dy );
			h2_kq_sqrt->Fill( i, sqrt(dq*dq + sq*sq) );
			h2_kr_sqrt->Fill( i, sqrt(dx*dx + dy*dy) );
		} // end-of-loop over FOOT layers.

		const mnd::geom::Line3D ld = RNTrackToLine3D( t );
		const mnd::geom::Point2D r0_down = ld.Eval( TFOOTHitCont::TARGET_Z );

		const double dx_dstr_ustr = r0_down.x - x0_upst;
		const double dy_dstr_ustr = r0_down.y - y0_upst;
		
		h2_ktxy->Fill(dx_dstr_ustr, dy_dstr_ustr);
		h1_kt->Fill(dx_dstr_ustr*dx_dstr_ustr + dy_dstr_ustr*dy_dstr_ustr);
		h1_kt_sqrt->Fill(sqrt(dx_dstr_ustr*dx_dstr_ustr + dy_dstr_ustr*dy_dstr_ustr));

		h1_diff_upstr_down-> Fill( ld.DistanceTo(lu) );
		h1_score->Fill(t.score);

		h1_sci21_cut2->Fill(sci21.E);
		h1_sci22_cut2->Fill(sci22.E);
		h1_sci31_cut2->Fill(sci31.E);
	}

	std::vector<std::string> resolutions {};
	
	/* Around each residue, also fit a teeny-weeny gauss-chan 🥺 👉👈
	 * Its' width is what we label as detector resolution. */
	TCanvas* c1d = new TCanvas("FitResidue1D", "Fit residues 1D", 2150, 1650);
	c1d->Divide(N_PAIRS, 3);
	for(size_t i=0; i<N_PAIRS; ++i) {
		c1d->cd(i+1);
		auto [fit_result_x, _x_] = resx[i]->DrawAndFit(2.0, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu X: " KBH_MAG "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_x[1], 1e3*fit_result_x[2]);
		c1d->cd(i+1 + N_PAIRS);
		auto [fit_result_y, _y_] = resy[i]->DrawAndFit(2.0, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu Y: " KBH_YEL "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_y[1], 1e3*fit_result_y[2]);

		c1d->cd(i+1 + 2*N_PAIRS);
		auto [fit_result_q, _q_] = resq[i]->DrawAndFit(2.0, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu Q: " KBH_CYN "%.3f, %.3f" KRNM " [Q]\n", i, fit_result_q[1], fit_result_q[2]);

		resolutions.emplace_back (
			Form("FOOT%zu (X,Y): (%s%.2f%s, %s%.2f%s) um",
				i, KRED, 1e3*fit_result_x[2], KNRM, KBLU, 1e3*fit_result_y[2], KNRM)
		);
	}
	TCanvas* cpos = new TCanvas("foot_xyq", "FOOT Positions/charge", 2150, 1650);
	cpos->Divide(N_PAIRS, 3);
	for(size_t i=0; i<N_PAIRS; ++i) {
		cpos->cd(i+1+ 0*N_PAIRS);
		h1_footx[i]->Draw();
		cpos->cd(i+1+ 1*N_PAIRS);
		h1_footy[i]->Draw();
		cpos->cd(i+1+ 2*N_PAIRS);
		h1q[i]->Draw();
	}
	
	TCanvas* c2d = new TCanvas("FitResidue2D", "Fit residues 2D", 2200, 1400);
	c2d->Divide(1, 3);
	c2d->cd(1); h2_resx->Draw("COLZ");
	c2d->cd(2); h2_resy->Draw("COLZ");
	c2d->cd(3); h2_resq->Draw("COLZ");

	[[maybe_unused]]
	TCanvas* cq2d = new TCanvas("Charge2D", "Charge of track measured by layers", 2000, 1200);
	h2q->Draw("COLZ");
	
	TCanvas* ckv = new TCanvas("CostValues2D", "Cost Values 2D", 2200, 1400);
	ckv->Divide(4, 2);
	ckv->cd(1); h2_kq->Draw("COLZ");
	ckv->cd(5); h2_kq_sqrt->Draw("COLZ");
	ckv->cd(2); h2_kr->Draw("COLZ");
	ckv->cd(6); h2_kr_sqrt->Draw("COLZ");
	ckv->cd(3); h1_kt->Draw();
	ckv->cd(7); h1_kt_sqrt->Draw();

	ckv->cd(4); h2_ktxy->Draw("COLZ");
	ckv->cd(8); h1_diff_upstr_down->Draw();

	TCanvas* cscore = new TCanvas("Score", "Score of the recognized tracks", 1200, 800);
	cscore->Divide(2,2);
	cscore->cd(1); h1_score->Draw();
	cscore->cd(2); new PLatex(0.08,
		"Coefficients: ",
		Form("Cr = %.1f", Cr),
		Form("Cq = %.1f", Cq),
		Form("Ct = %.1f", Ct),
		Form("max cost: %.1f", max_cost),
		Form("max cost_f: %.1f", max_cost_f)
	);
	cscore->cd(3); h1_has_upstream->Draw();
	cscore->cd(4); new PLatex(0.08,
		"Position resolution: ",
		mnd::as_span(resolutions)
	);

	TCanvas* cR = new TCanvas("RScore", "Individual R-score component", 2150, 1400);
	TCanvas* cQ = new TCanvas("QScore", "Individual Q-score component", 2150, 1400);
	TCanvas* cT = new TCanvas("TScore", "Individual T-score component", 2150, 1400);
#define DECORATE_HISTO(hshort, col) \
	for(u32 n=0; n<N_PAIRS; ++n) { \
		TH1I* h = h1_diff_##hshort.at(n); \
		h->SetFillColor(col); h->SetLineColor(kBlack); h->SetLineWidth(1);  \
		h = h1_acc_##hshort.at(n); \
		h->SetFillColor(col); h->SetLineColor(kBlack); h->SetLineWidth(1);  \
	}
	DECORATE_HISTO(r, kRed-7);
	DECORATE_HISTO(q, kCyan-4);
	DECORATE_HISTO(t, kYellow-7);
	cR->Divide(4,2);
	cQ->Divide(4,2);
	cT->Divide(4,2);
	static_assert(N_PAIRS == 4);
	for(u32 n=0; n<N_PAIRS; ++n) {
		cR->cd(1 + n); h1_diff_r.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
        cR->cd(5 + n); h1_acc_r.at( n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
		cQ->cd(1 + n); h1_diff_q.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
        cQ->cd(5 + n); h1_acc_q.at( n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
		cT->cd(1 + n); h1_diff_t.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
        cT->cd(5 + n); h1_acc_t.at( n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
	}

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2150, 1400);
	cs->Divide(3,3);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	cs->cd(7); h1_sci21_cut2->Draw();
	cs->cd(8); h1_sci22_cut2->Draw();
	cs->cd(9); h1_sci31_cut2->Draw();

	double qvg = cut_q.size() > 0
	? (std::accumulate(cut_q.begin(), cut_q.end(), 0.0,
		[](double sum, const auto& value) -> double {
			return sum + (value.second[0] + value.second[1]) / 2;
		}) / cut_q.size() )
	: 0.0;

	canvas::save_all<canvas::Exe>(save, {
		mnd::file_stem(fileName),
		std::string{"Q:"} + mnd::ftos(qvg)
	});

	WARN("End-of-main");
	rootApp.Run();

#else
    WARN("\n" BOLD ".. MND_FOOTTRACK_DEBUG not compiled in. This program does nothing.\n" KNRM);
#endif
    return 0;
}
