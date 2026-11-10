#include "util/CLI.h"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/FFT.h"
#include "util/Option.hxx"
#include "util/TextUtils.h"

#include "ExperimentAssumptions.hh"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

using mnd::Option;
using DoFit = Option<std::vector<i32>>;

constexpr f64 DELTA_EDGE_DEFAULT = 0.4999;
constexpr i32 Z0 = (i32)mnd::assume::primary.Z;

/* Because cling issues the WEIRDEST compiler error,..
 * I cannot just simply return `FOOTDeltaParam` instance. It tries to compile the class
 * from the inputs but loses it on template spec in boost preproc library LOL. PEGI 18. */
class FOOTDeltaParam;
FOOTDeltaParam* GetDeltaParams(const TH2D& , u32);

int main(int argc, char* argv[]) {
	CLI::App app{"Perform the delta-correction of a specific FOOT detector. "
		"NOTE: to be valid, must be done *after* the gain-matching, and requires proper gain-matched values. "};

	std::vector<std::string> fileNames {};
	i32 ifoot = 0;
	A3 foot_binning = {80, 1600,4600};
	A3 delta_binning = {80, -DELTA_EDGE_DEFAULT, DELTA_EDGE_DEFAULT};
	i32 Q_target = 6;
	
	auto sratio = 0.9_f64;
	auto niter  = 2_u32;
	auto harmonics = 8_i32;
	Verbosity v = Verbosity::INFO;

	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	A3 sci_binning = {500, 300, 4000};
	
	std::vector<canvas::Extension> save{};
	bool plot_gain = false;
	DoFit do_fit {};

	u32 nthreads = 1;

	add_logged_option(app, "-f,--file", fileNames, "Pass a file name.")
		->check(CLI::ReadPermissions)
        ->delimiter(',');
	add_logged_option(app, "-i,--foot-id", ifoot,
		"Select which FOOT detector.")
		->check(CLI::Range(0, N_FOOT_DETECTORS-1));
	add_logged_option(app, "-q,--charge", Q_target,
		"Select which target charge we're roughly gating upon.")
		->check(CLI::Range(1, Z0));
	
	add_logged_option(app, "-b,--foot-binning", foot_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "-d,--delta-binning", delta_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "-l,--fit", do_fit, "Select which ASICs to fit.")
		->delimiter(',')
		->type_name("INT[,INT...]")
		->check(CLI::RangeOrEmpty(0,10));
	add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--harmonics", harmonics, "Select total number of harmonics for the FFT.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--verbosity", v, "Select verbosity level of the FFT fitter.");
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "--sci-binning",sci_binning, "Scintillator QDC binning")
		->delimiter(',');
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');

	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
	)->check(CLI::PositiveNumber);
	add_logged_flag(app, "-g,--plot-gain,!--no-plot-gain", plot_gain, "Do plot gain");
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.size() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, 0,
			i * 64, 1000);
		line->SetLineColor(kBlack);
		line->SetLineStyle(2);
		line->SetLineWidth(2);
		vlines.push_back( line );
	}
#define DRAW_VLINES(name) \
	{ \
		TH2D* h2_ = name; \
		for(auto* l0 : vlines) { \
			TLine* l = dynamic_cast<TLine*>(l0->Clone()); \
			l->SetY1(h2_->GetYaxis()->GetXmin()); \
			l->SetY2(h2_->GetYaxis()->GetXmax()); \
			l->Draw("SAME"); \
		} \
		gPad->SetGridx(false); gPad->Update(); \
	}

	FOOTParam* p;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileNames.front().c_str(), "READ");
		get_obj(f, p, Form("FOOT%d_setup", ifoot));
	}

	const auto& gain = p->gain;
	
	/* Check if in gain we have the entry for Z=Q_target. */
	gain.IsSaneThrowable();

	const double CA = gain.GetNominalValue(Q_target);
	if(!std::isfinite(CA))
		ERROR("Z=%d nominal not found in the setup file. Cannot proceed.\n", Q_target);

	ROOT::EnableThreadSafety();
	ROOT::DisableImplicitMT();
	TApplication rootApp("app", 0, 0);

	constexpr A3 x_binning = {160, 0, RNFOOTCal::N_STRIPS};
	TH1P h1_delta = TH1P(Form("((h1_foot%d))Delta@FOOT%d", ifoot, ifoot), kGreen+1, delta_binning);
	TH2P sum_energy_vs_x = TH2P(Form("((h2_footraw%d))Cluster sum [ADC]:Strip num@FOOT%d gain matched", ifoot, ifoot),
		x_binning, foot_binning);
	TH2P sum_energy_vs_delta = TH2P(Form("((h2_footraw%d))Cluster sum [ADC]:Delta@FOOT%d gain matched", ifoot, ifoot),
		delta_binning, foot_binning);
	TH2P corr_energy_vs_delta = TH2P(Form("((h2_footcorr%d))E1/%.1f - 1:Delta@FOOT%d", ifoot, CA, ifoot),
		delta_binning,
		foot_binning[0], foot_binning[1]/CA - 0.9, foot_binning[2]/CA - 0.9);

	auto h1_m0 = TH1P(Form("((h_0)) FOOT%d dE [ADC units]@All cluster sizes", ifoot), kYellow - 7, foot_binning);

	auto h1_sci21 = TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, sci_binning);
	auto h1_sci22 = TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, sci_binning);
	auto h1_sci31 = TH1P("SCI31 QDC mean [QDC units]", 0x009B2F_c, sci_binning);
	auto h1_sci21_cut = TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, sci_binning);
	auto h1_sci22_cut = TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, sci_binning);
	auto h1_sci31_cut = TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, sci_binning);
	
    show_console_cursor(false);

	DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
        auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fname);

        const size_t nentries = ntuple->GetNEntries();

		auto _bar_ptr = std::make_unique<ProgressBar>(
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg(
					"%zu/%zu: %'zu (%s)", i+1, fileNames.size(), nentries, fname.c_str()
				)},
                option::ForegroundColor{ indicators::next_col() },
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];

        for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);

			const auto& sci21 = frs->sci[0];
            const auto& sci22 = frs->sci[1];
            const auto& sci31 = frs->sci[2];

            if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
            if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
            if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);

            if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
            if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
            if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);
			
			for(const auto& cl : foot->fCl) {
				if(cl.fCM == 1) continue;

				double delta = cl.Delta();
				double e = cl.fCE;
				double x = cl.fCX;
				
				e *= gain.CorrectionFactor(x, e);

				/* Explicitly skip here so we don't deal with over/underflow bins */
				if(!mnd::IsInside(e, foot_binning)) continue;

				sum_energy_vs_delta -> Fill(delta, e);
				sum_energy_vs_x -> Fill(x, e);
				h1_delta -> Fill(delta);
				h1_m0->Fill(e);
			}
		}
		if(!bar.is_completed())
			bar.mark_as_completed();
    });

	/* Perform the Step (1) of the delta correction.
	 * Of course, this part isn't threadsafe, so yeet it out of the paralle_process. */
	auto* d = GetDeltaParams(sum_energy_vs_delta, niter);
	WARN("\nRegoing over the file ... Woohoooo \n");
	
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fname);

        const size_t nentries = ntuple->GetNEntries();

		auto _bar_ptr = std::make_unique<ProgressBar>(
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg(
					"%zu/%zu: %'zu (%s)", i+1, fileNames.size(), nentries, fname.c_str()
				)},
                option::ForegroundColor{ indicators::next_col() },
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			
			const auto& sci21 = frs->sci[0];
			const auto& sci22 = frs->sci[1];
			const auto& sci31 = frs->sci[2];
			
			if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);

			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);

			for(const auto& cl : foot->fCl) {
				if(cl.fCM == 1) continue;

				double delta = cl.Delta();
				double e = cl.fCE;
				double x = cl.fCX;
				e *= gain.CorrectionFactor(x, e);

				/* Normalize with the triangular delta-correction */
				e /= d->CorrectionBasic(delta);
				
				corr_energy_vs_delta -> Fill(delta, e / CA - 1);
			}
		}
		if(!bar.is_completed())
			bar.mark_as_completed();
	});

    show_console_cursor(true);

	/* If: f(x) ≈ a0 + sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * Then e2 := e1 - sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * to fit everything onto the DC line. */
	auto [fft, gr, g] = DoFFTW<fit_info::GAUSS_MAX>(corr_energy_vs_delta, -0.5, 0.5, sratio, harmonics, v);
	WARN("\n[[FOOT%d]] FFT Params:\n", p->de10_index_);
	printf("\"s\": %.5f,\n\"f\": %.5f,\n", d->s, d->f);
	printf("\"f2\": {\n\t\"n\": %d,", static_cast<i32>(delta_binning[0]));
	printf("\n\t\"c\": [");
	for(int i=0; i<harmonics; ++i) std::cout << fft.coeff[i] << ", ";
	std::cout << fft.coeff[harmonics] << "]";
	printf("\n}\n");

	TCanvas *c = new TCanvas("FOOT", Form("Delta%d", ifoot), 2000, 1400);
	c->Divide(3,2);
	c->cd(1);
	h1_delta->Draw();
	c->cd(2);
	// h2_m_vs_delta->Draw("COLZ");
	c->cd(3); gPad->SetLogz();
	sum_energy_vs_delta->Draw("COLZ");
	c->cd(4);  gPad->SetLogz();
	corr_energy_vs_delta->Draw("COLZ");
	g->Draw("L SAME"); gPad->SetLogz();
	// gr->Draw("P SAME");
	c->cd(5);
	sum_energy_vs_x->Draw("COLZ");
	c->cd(6);
	h1_m0->Draw();

	TCanvas* cs = new TCanvas("SCI1D", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* cs2 = new TCanvas("SCI2D", "SCI-Corr", 2200, 1400);
	cs2->Divide(2,2);
	// cs2->cd(1); h2_sci2v2->Draw("COLZ");
	// cs2->cd(2); h2_sci3v2->Draw("COLZ");
	// cs2->cd(3); h2_sci2v2_cut->Draw("COLZ");
	// cs2->cd(4); h2_sci3v2_cut->Draw("COLZ");

	if(plot_gain) {
		TCanvas* gn = new TCanvas("gain", "Gain factor", 2000, 1300);
		gn->Divide(2,2);
		TH2D* hgain = gain.GetHisto();
		hgain->SetTitle(Form("FOOT%d Gain Parameter", ifoot));
		gn->cd(1); hgain->Draw("COLZ"); DRAW_VLINES(hgain)
		gn->cd(2); hgain->Draw("SURF1");
		auto [p4, graph4] = gain.GetGraph(4*64 + 32);
		gn->cd(3); graph4->Draw("AL"); p4->Draw("P SAME");
		new PLatex(0.08, "Middle of ASIC[4]");
		auto [p5, graph5] = gain.GetGraph(5*64 + 32);
		gn->cd(4); graph5->Draw("AL"); p5->Draw("P SAME");
		new PLatex(0.08, "Middle of ASIC[5]");
	}

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOT%d", ifoot),
		Form("Z_%d",  Q_target),
		mnd::rootfile_names_concatenated(fileNames)
	});
	WARN("End-of-main\n");
	rootApp.Run(); return 0;
}

FOOTDeltaParam* GetDeltaParams(const TH2D& e_vs_delta, u32 niter) {
	constexpr static std::array<double, 2> SIDE_PEAK = {0.2, 0.4};
	double sp, sm, A0, Ap, Am;
	int bmax_0, bmax_p, bmax_m;
	auto* result = new FOOTDeltaParam{};

	auto h1_delta = std::unique_ptr<TH1D>( e_vs_delta.ProjectionX("px_") );
	h1_delta->SetDirectory(nullptr);
	
	TAxis* ax = h1_delta->GetXaxis();

	ax->SetRangeUser( -SIDE_PEAK[0], SIDE_PEAK[0] );
	bmax_0 = h1_delta->GetMaximumBin();
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( SIDE_PEAK[0], SIDE_PEAK[1] );
	bmax_p = h1_delta->GetMaximumBin();
	sp = ax->GetBinCenter(bmax_p);
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( -SIDE_PEAK[1], -SIDE_PEAK[0] );
	bmax_m = h1_delta->GetMaximumBin();
	sm = ax->GetBinCenter(bmax_m);
	ax->SetRange(0,0); // unzoom
	double s = (sp - sm)/2;
	result->s = s;

	/* Next, in the TH2D just try to project the fews bins around sp,sm. */

	auto h1_e_0 = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_0", bmax_0-1, bmax_0+1) );
	auto h1_e_p = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_p", bmax_p-1, bmax_p+1) );
	auto h1_e_m = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_m", bmax_m-1, bmax_m+1) );
	h1_e_0->SetDirectory(nullptr);
	h1_e_p->SetDirectory(nullptr);
	h1_e_m->SetDirectory(nullptr);

	double side_ratio = 1;
	A0 = GaussFitMax( h1_e_0.get(), side_ratio, niter, Verbosity::CHATTY ).first.at(1);
	Ap = GaussFitMax( h1_e_p.get(), side_ratio, niter, Verbosity::CHATTY ).first.at(1);
	Am = GaussFitMax( h1_e_m.get(), side_ratio, niter, Verbosity::CHATTY ).first.at(1);

	double f = (Ap + Am) / (2 * A0);
	result->f = f;

	WARN("\nA0 = %.2f, A+ = %.2f, A- = %.2f\n"
	       "s0 = %.2f, s+ = %.2f, s- = %.2f\n"
		   ">> f = %.4f; s = %.4f <<\n"
		   "b0 = %d, b+ = %d, b- = %d\n",
		   A0,Ap,Am, 0.0,sp,sm, f,s, bmax_0,bmax_p,bmax_m);

	return result;
}
