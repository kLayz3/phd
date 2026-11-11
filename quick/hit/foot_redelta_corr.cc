/* This include is just to have LSP support when the define is present at compile-time
 * but not during editing. */
//#define MND_DEBUG_ENABLED

#include "util/CLI.h"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/FFT.h"
#include "util/Option.hxx"
#include "util/TextUtils.h"
#include "util/GaussFitMax.hxx"

#include "ExperimentAssumptions.hh"

#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"
#include "util/Verbosity.hxx"

using mnd::Option;
using mnd::None;

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

constexpr auto N_PAIRS = TFOOTHitCont::N_PAIRS;

constexpr f64 DELTA_EDGE_DEFAULT = 0.49999;
constexpr i32 Z0 = (i32)mnd::assume::primary.Z;

class FOOTDeltaParam;
FOOTDeltaParam* GetDeltaParams(const TH2D& , u32 , Verbosity );

/* This procedure basically redoes the FOOT delta correction but instead of using cal level,
 * goes to hit level to be also able to gate on straight tracks. */
int main(int argc, char* argv[]) {
	ROOT::EnableThreadSafety();

	TClass* _cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!_cl or !_cl->GetDataMember("_u"))
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	CLI::App app{"Perform (re)delta-correction of a specific FOOT detector, "
	    "based on already preliminary results from the Kalman filter. Basically, here we can gate explicitly "
        "on specific 'Z' tracks (can be wrong z tho)."
		"NOTE: to be valid, must be done *after* the gain-matching, and requires proper gain-matched values."
	};

	std::vector<std::string> fileNames {};
	i32 ifoot = 0;
	A3 foot_binning = {80, 1600,4600};
	A3 delta_binning = {80, -DELTA_EDGE_DEFAULT, DELTA_EDGE_DEFAULT};
	
	i32 Q_target = 6;
	A2 Q_width = {-0.5, 0.5};
	bool require_all = false;
	
	A3 foot_q_binning = {80,0,7.5};
	
	auto sratio = 0.9_f64;
	auto niter  = 2_u32;
	auto harmonics = 8_i32;
	Verbosity v = Verbosity::INFO;

	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	A3 sci_binning = {500, 300, 4000};
	
	std::vector<canvas::Extension> save{};
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
	add_logged_option(app, "-w,--qwidth", Q_width, "Select width of the target charge to gate upon.")
		->delimiter(',');
	add_logged_flag(app, "-a,--all", require_all,
		mnd::msg("Require all %u layers to participate in the selected tracks.", N_PAIRS));
	add_logged_option(app, "-b,--foot-binning", foot_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "-d,--delta-binning", delta_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "--foot-q-binning",foot_q_binning, "FOOT track charge (Q) binning.")
		->delimiter(',');
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
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.size() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

#ifdef MND_FOOTTRACK_DEBUG
	const A2 Q_interval = { Q_target + Q_width[0], Q_target + Q_width[1] };

	const int foot_pair_id = ifoot / 2;
	std::unique_ptr<std::array<FOOTParam, 2>> setup;
	{
		auto f = std::make_unique<TFile>(fileNames.front().c_str(), "READ");
		std::array<FOOTParam, 2>* parameters = nullptr;
		get_obj(f, parameters, Form("FOOT_%d_setup", foot_pair_id));
		setup.reset(parameters);
	}
	const auto parameter = std::find_if(setup->begin(), setup->end(), [ifoot](const FOOTParam& p) {
		return p.N == ifoot;
	});
	if(parameter == setup->end())
		ERROR("FOOT%d not found in saved pair %d setup.\n", ifoot, foot_pair_id);
	
	FOOTParam* const foot_param = &*parameter;

	const auto& gain = foot_param->gain;
	const auto& delta_p = foot_param->de; // currently the delta-corr parametrisation.
	const Orientation o = foot_param->GetOrientation();

	/* Check if in gain we have the entry for Z=Q_target. */
	gain.IsSaneThrowable();

	if(harmonics > static_cast<i32>(delta_binning[0]) / 2)
		ERROR("Cannot select %d harmonics, which is >%d  (binning/2).\n",
			harmonics, static_cast<i32>(delta_binning[0]) / 2);

	const double CA = gain.GetNominalValue(Q_target);
	if(!std::isfinite(CA))
		ERROR("Z=%d nominal not found in the setup file. Cannot proceed.\n", Q_target);

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

	auto h1_foot_q  = TH1P{"FOOT Track Q [e]", 0xFF7C0A_c, foot_q_binning};
	auto h1_foot_qc = TH1P{"FOOT Track Q [e]", 0xFFAD64_c, foot_q_binning};

	show_console_cursor(false);

	DynamicProgress<ProgressBar> progress1;
	mnd::parallel_process(fileNames, nthreads, [=, &progress1](size_t i, auto fname) mutable {
        auto model = RNTupleModel::Create();
        auto foot = model->MakeField<RNFOOTHit>("FOOT");
        auto frs = model->MakeField<RNFRSHit>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);

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
		auto idx = progress1.push_back(std::move(_bar_ptr));
		auto& bar = progress1[idx];

        for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
            if(mnd::PrintProgress(bar, entryId, nentries, 1000))
				progress1.print_progress();

			const auto& sci21 = frs->cal.sci[0];
            const auto& sci22 = frs->cal.sci[1];
            const auto& sci31 = frs->cal.sci[2];

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
			
			for(const auto& t : foot->track) {
				h1_foot_q.Fill(t.Q);

                if(!mnd::IsInside(t.Q, Q_interval)) continue;
				if(require_all && t.n != N_PAIRS) continue;
				
				/* All of the cal fields are behind some kind of a `T: std::array<T, N_FOOT_PAIRS>` */
                const double delta = (o == Orientation::X)
                    ? t._delta_x[foot_pair_id]
                    : t._delta_y[foot_pair_id];

                /* nil value is always mapped to a quiet nan */
                if(!std::isfinite(delta)) continue;
				const double m = (o == Orientation::X)
					? t._u[foot_pair_id]
					: t._v[foot_pair_id];

				/* m is the 'bare' measurement. */
				const double cx = foot_param->BareToStripnum(m);
                const double e0 = (o == Orientation::X)
                    ? t._e0_x[foot_pair_id]
                    : t._e0_y[foot_pair_id];
				
				const double gm = foot_param->gain.CorrectionFactor(cx, e0);
				const double e = e0 * gm;
				
				sum_energy_vs_delta -> Fill(delta, e);
				sum_energy_vs_x -> Fill(cx, e);
				h1_delta -> Fill(delta);
				h1_m0->Fill(e);
				
				h1_foot_qc.Fill(t.Q);
			}
		}
		if(!bar.is_completed())
			bar.mark_as_completed();
    });

	/* Perform the Step (1) of the delta correction.
	 * Of course, this part isn't threadsafe, so yeet it out of the paralle_process. */
	const auto* d = GetDeltaParams(sum_energy_vs_delta, niter, v);
	WARN("\nRegoing over the file ... Woohoooo \n");
	
	DynamicProgress<ProgressBar> progress2;
	mnd::parallel_process(fileNames, nthreads, [=, &progress2](size_t i, auto fname) mutable {
        auto model = RNTupleModel::Create();
        auto foot = model->MakeField<RNFOOTHit>("FOOT");
        auto frs = model->MakeField<RNFRSHit>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);

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
		auto idx = progress2.push_back(std::move(_bar_ptr));
		auto& bar = progress2[idx];

        for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
            if(mnd::PrintProgress(bar, entryId, nentries, 1000))
				progress2.print_progress();

			const auto& sci21 = frs->cal.sci[0];
            const auto& sci22 = frs->cal.sci[1];
            const auto& sci31 = frs->cal.sci[2];

            if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
            if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
            if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;
			
			for(const auto& t : foot->track) {
                if(!mnd::IsInside(t.Q, Q_interval)) continue;
				if(require_all && t.n != N_PAIRS) continue;
				
				/* All of the cal fields are behind some kind of a `T: std::array<T, N_FOOT_PAIRS>` */
                const double delta = (o == Orientation::X)
                    ? t._delta_x[foot_pair_id]
                    : t._delta_y[foot_pair_id];

                /* nil value is always mapped to a quiet nan */
                if(!std::isfinite(delta)) continue;
				const double m = (o == Orientation::X)
					? t._u[foot_pair_id]
					: t._v[foot_pair_id];

				/* m is the 'bare' measurement. */
				const double cx = foot_param->BareToStripnum(m);
                const double e0 = (o == Orientation::X)
                    ? t._e0_x[foot_pair_id]
                    : t._e0_y[foot_pair_id];
				
				const double gm = foot_param->gain.CorrectionFactor(cx, e0);
				double e = e0 * gm;
				
				/* Normalize with the triangular delta-correction */
				e /= d->CorrectionBasic(delta);
				
				corr_energy_vs_delta->Fill(delta, e / CA - 1);
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
	WARN("\n[[FOOT%d]] FFT Params:\n", foot_param->de10_index_);
	printf("\"s\": %.6f,\n\"f\": %.6f,\n", d->s, d->f);
	printf("\"f2\": {\n\t\"n\": %d,", static_cast<i32>(delta_binning[0]));
	printf("\n\t\"c\": [");
	for(int i=0; i<harmonics; ++i) std::cout << fft.coeff[i] << ", ";
	std::cout << fft.coeff[harmonics] << "]";
	printf("\n}\n");

	TCanvas *cmain = new TCanvas("Delta", Form("Delta-correction: FOOT%d", ifoot), 2000, 1400);
	cmain->cd(); corr_energy_vs_delta->Draw("COLZ"); gPad->SetLogz();
	g->SetLineColor(kRed + 1);
	g->SetLineStyle(kSolid);
	g->SetLineWidth(4);
	g->Draw("L SAME");
	gr->SetMarkerColor(kBlack);
	gr->SetLineColor(kBlack);
	gr->Draw("P SAME");
	TGraph* fft_old = delta_p.f2.GetGraph();
	fft_old->SetLineColor(kBlue + 1);
	fft_old->SetLineStyle(kDashed);
	fft_old->SetLineWidth(4);
	fft_old->Draw("L SAME");
	auto* fft_legend = new TLegend(0.12, 0.75, 0.48, 0.9);
	fft_legend->AddEntry(gr, "Gaussian slice means", "pe");
	fft_legend->AddEntry(g, "New Fourier correction", "l");
	fft_legend->AddEntry(fft_old, "Old Fourier correction (stored setup)", "l");
	fft_legend->Draw();

	TCanvas *c = new TCanvas("FOOT", Form("Delta%d", ifoot), 2000, 1400);
	c->Divide(2,2);
	c->cd(1); h1_delta->Draw();
	c->cd(2); sum_energy_vs_delta->Draw("COLZ"); gPad->SetLogz();
	TGraph* graph_basic = d->GetGraphBasic(CA);
	TGraph* graph_old = delta_p.GetGraphBasic(CA);

	graph_basic->SetLineColor(kRed + 1);
	graph_basic->SetLineStyle(kSolid);
	graph_basic->SetLineWidth(4);
	graph_basic->Draw("L SAME");
	graph_old->SetLineColor(kBlue + 1);
	graph_old->SetLineStyle(kDashed);
	graph_old->SetLineWidth(4);
	graph_old->Draw("L SAME");
	auto* basic_legend = new TLegend(0.12, 0.78, 0.56, 0.9);
	basic_legend->AddEntry(graph_basic, "New triangular correction", "l");
	basic_legend->AddEntry(graph_old, "Old triangular correction (stored setup)", "l");
	basic_legend->Draw();

	c->cd(3); sum_energy_vs_x->Draw("COLZ");
	c->cd(4); h1_m0->Draw();

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas *cf = new TCanvas("ft", "FOOT Track Selection", 1400,800);
	cf->Divide(2,1);
	cf->cd(1); h1_foot_q.Draw();
	cf->cd(2); h1_foot_qc.Draw();

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOT%d", ifoot),
		Form("Z_%d",  Q_target),
		mnd::rootfile_names_concatenated(fileNames)
	});
	WARN("End-of-main\n");
	rootApp.Run();

#else
    WARN("\n" BOLD ".. MND_FOOTTRACK_DEBUG not compiled in. This program does nothing.\n" KNRM);
#endif
    return 0;

}

FOOTDeltaParam* GetDeltaParams(const TH2D& e_vs_delta, u32 niter, Verbosity v) {
	/* Interval where the 'side' peak is expected. */
	constexpr A2 SIDE_PEAK = {0.2, 0.4};
	constexpr f64 SRATIO = 1.5;

	double s0, sp, sm;
	double A0, Ap, Am;

	auto* result = new FOOTDeltaParam{};

	auto h1_delta = std::unique_ptr<TH1D>( e_vs_delta.ProjectionX("px_") );
	h1_delta->SetDirectory(nullptr);
	
	TAxis* ax = h1_delta->GetXaxis();

	ax->SetRangeUser( -SIDE_PEAK[0], SIDE_PEAK[0] );
	auto results_0 = GaussFitMax(h1_delta.get(), SRATIO, 5u, Verbosity::SILENT, true);
	if(!VerifyGaussFitMax(results_0) ) {
		WARN(">>> GaussFitMax invalid results around 0 peak; falling back to simple peak position.\n");
		s0 = h1_delta->GetBinCenter( h1_delta->GetMaximumBin() );
	} else {
		s0 = results_0.first[1]; // μ
	}
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( SIDE_PEAK[0], SIDE_PEAK[1] );
	auto results_p = GaussFitMax(h1_delta.get(), SRATIO, 5u, Verbosity::SILENT, true);
	if(!VerifyGaussFitMax(results_p) ) {
		WARN(">>> GaussFitMax invalid results around + (plus) peak; falling back to simple peak position.\n");
		sp = h1_delta->GetBinCenter( h1_delta->GetMaximumBin() );
	} else {
		sp = results_p.first[1]; // μ
	}
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( -SIDE_PEAK[1], -SIDE_PEAK[0] );
	auto results_m = GaussFitMax(h1_delta.get(), SRATIO, 5u, Verbosity::SILENT, true);
	if(!VerifyGaussFitMax(results_m) ) {
		WARN(">>> GaussFitMax invalid results around - (minus) peak; falling back to simple peak position.\n");
		sm = h1_delta->GetBinCenter( h1_delta->GetMaximumBin() );
	} else {
		sm = results_m.first[1]; // μ
	}
	ax->SetRange(0,0); // unzoom
	
	double s = (sp - sm)/2;
	result->s = s;

	/* Next, in the TH2D just try to project the fews bins around sp,sm. */
	const int bmax_0 = h1_delta->FindBin(s0);
	const int bmax_p = h1_delta->FindBin(sp);
	const int bmax_m = h1_delta->FindBin(sm);

	auto h1_e_0 = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_0", bmax_0-1, bmax_0+1) );
	auto h1_e_p = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_p", bmax_p-1, bmax_p+1) );
	auto h1_e_m = std::unique_ptr<TH1D>( e_vs_delta.ProjectionY("py_m", bmax_m-1, bmax_m+1) );
	h1_e_0->SetDirectory(nullptr);
	h1_e_p->SetDirectory(nullptr);
	h1_e_m->SetDirectory(nullptr);

	double side_ratio = 1;
	const auto fit0 = GaussFitMax( h1_e_0.get(), side_ratio, niter, v);
	const auto fitp = GaussFitMax( h1_e_p.get(), side_ratio, niter, v);
	const auto fitm = GaussFitMax( h1_e_m.get(), side_ratio, niter, v);

	if(!VerifyGaussFitMax(fit0))
		ERROR("Central fit slice has non-convergent Gauss fit.");
	if(!VerifyGaussFitMax(fitp))
		ERROR("Positive-side fit slice has non-convergent Gauss fit.");
	if(!VerifyGaussFitMax(fitm))
		ERROR("Negative-side fit slice has non-convergent Gauss fit.");

	A0 = fit0.first.at(1);
	Ap = fitp.first.at(1);
	Am = fitm.first.at(1);

	double f = (Ap + Am) / (2 * A0);
	result->f = f;

	WARN("\nA0 = %.4f, A+ = %.4f, A- = %.4f\n"
	       "s0 = %.4f, s+ = %.4f, s- = %.4f\n"
		   ">> f = %.6f; s = %.6f <<\n"
		   "b0 = %d, b+ = %d, b- = %d\n",
		   A0,Ap,Am, s0,sp,sm, f,s, bmax_0,bmax_p,bmax_m);

	return result;
}
