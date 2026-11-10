/* This include is just to have LSP support when the define is present at compile-time
 * but not during editing. */
//#define MND_DEBUG_ENABLED

#include "util/CLI.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/FitSpline.h"
#include "util/TextUtils.h"

#include "TStyle.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"

enum class Take { gauss, profile, gauss_fit_only };

using mnd::Option;
using mnd::None;

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

constexpr auto N_PAIRS = TFOOTHitCont::N_PAIRS;

/* This procedure basically redoes the FOOT gain match but instead of using cal level,
 * goes to hit level to be also able to gate on straight tracks. */
int main(int argc, char* argv[]) {
	ROOT::EnableThreadSafety();
	TClass* _cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!_cl or !_cl->GetDataMember("_u"))
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	CLI::App app{"Perform (re)gain matching of a specific FOOT detector, "
	    "based on already preliminary results from the Kalman filter. Basically, here we can gate explicitly "
        "on specific 'Z' tracks (can be wrong z tho)"};
	
    std::vector<std::string> fileNames;
	int ifoot = 0;
	int bins_per_asic = 64;
	Option<std::vector<i32>> do_fit = None;
	double sratio = 0.9;
	u32 niter = 2;
	A3 foot_binning = {1000, 4, 4000};
	A3 delta_binning = {150, -0.499, 0.499};
	A3 foot_q_binning = {80,0,7.5};
	A3 sci_binning = {500, 300, 4000};
	double delta_cut = 0.05;
	int Q_target = 6;
	A2 Q_width = {-0.5,0.5};
	bool require_all = false;
	size_t poly_deg = 4;
	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	std::vector<canvas::Extension> save = {};

	u32 nthreads = 1;

	bool add_asic_indicator = false;
	Take take = Take::gauss_fit_only;
	Option<std::vector<i32>> show_old = None;

	add_logged_option_seq(app, "-f,--file", fileNames, "Pass one or more file name.")
		->check(CLI::ReadPermissions)
        ->delimiter(',');
	add_logged_option(app, "-i,--foot-id", ifoot,
		"Select which FOOT detector.")
		->check(CLI::Range(0, N_FOOT_DETECTORS-1));
	add_logged_option(app, "-q,--charge", Q_target,
		"Select which target charge to gate upon.")
		->check(CLI::Range(1,6));
	add_logged_option(app, "-w,--qwidth", Q_width,
		"Select width of the target charge to gate upon.")
		->delimiter(',');
	add_logged_flag(app, "-a,--all", require_all,
		mnd::msg("Require all %u layers to participate in the selected tracks.", N_PAIRS));
	add_logged_option(app, "-p,--poly", poly_deg,
		"Select polynomial degree to fit the specified ASIC's.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-d,--delta", delta_cut,
		"Select the delta cut which will be applied to gate on very central hits.")
		->check(CLI::Range(0.01, 0.2));
	add_logged_option(app, "--foot-binning",foot_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "--delta-binning",delta_binning, "δ- value binning")
		->delimiter(',');
	add_logged_option(app, "--foot-q-binning",foot_q_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "-b,--bins",bins_per_asic, "Binning Per ASIC")
		->check([](const std::string& match) -> std::string {
			int n = -1;
			if(!mnd::parse(match, n) or n <= 0 or (TFOOTCalCont::N_STRIPS % n != 0))
				return mnd::msg("Must be integer parsable and evenly divisible by %d", TFOOTCalCont::N_STRIPS);
			return "";
		});
	add_logged_option(app, "-l,--fit", do_fit, "Select which ASICs to fit.")
		->delimiter(',');
	add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "--sci-binning",sci_binning, "Scintillator QDC binning")
		->delimiter(',');
	add_logged_option(app, "-t,--take", take,
		"Which type of projection fit to take. `gauss_fit_only` will only use gaussian spline for well sampled data.");
	add_logged_option(app, "-s,--show-old", show_old,
		"Sequence of Z charges (, sep) to overlay their current gain match curve on the canvas.")
		->delimiter(',')
		->type_name("INT[,INT...]");
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
	)->check(CLI::PositiveNumber);
	add_logged_flag(app, "--add-asic-indicator", add_asic_indicator,
		"add an \"asic\": NUM indicator to the output JSON. For easier navigation.");

	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.size() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

#ifdef MND_FOOTTRACK_DEBUG
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
	
	const Orientation o = foot_param->GetOrientation();
	A2 delta_interval { -delta_cut, delta_cut };
	const A2 Q_interval = { Q_target + Q_width[0], Q_target + Q_width[1] };
	
	TApplication rootApp("app", 0, 0);

	auto h1_sci21     = TH1P{"SCI21 QDC mean [QDC units]", 0xCB00CB_c, sci_binning};
	auto h1_sci22     = TH1P{"SCI22 QDC mean [QDC units]", 0x0070DD_c, sci_binning};
	auto h1_sci31     = TH1P{"SCI31 QDC mean [QDC units]", 0x009B2F_c, sci_binning};
	auto h1_sci21_cut = TH1P{"((h1_cut)) SCI21 QDC mean [QDC units]@With cut and FOOT selection", 0x890389_c, sci_binning};
	auto h1_sci22_cut = TH1P{"((h1_cut)) SCI22 QDC mean [QDC units]@With cut and FOOT selection", 0x6180FD_c, sci_binning};
	auto h1_sci31_cut = TH1P{"((h1_cut)) SCI31 QDC mean [QDC units]@With cut and FOOT selection", 0x7DE69D_c, sci_binning};

	auto h1_foot_q  = TH1P{"FOOT Track Q [e]", 0xFF7C0A_c, foot_q_binning};
	auto h1_foot_qc = TH1P{"FOOT Track Q [e]", 0xFFAD64_c, foot_q_binning};

	auto h1_delta = TH1P{"((h1_d0))Delta [from -0.5, 0.5]", kGreen-2, delta_binning};

	auto h2_e0_vs_x = TH2P {
		Form("((h2_r))Cluster raw energy [ADC]:Strip number [0..640]@FOOT%d E raw, requested Q=%d", ifoot, Q_target),
		bins_per_asic*TFOOTCalCont::N_ASIC, 0, RNFOOTCal::N_STRIPS, foot_binning
	};

	auto h2_egm_vs_d = TH2P {
		Form("((h2_d))Cluster energy gain-matched [ADC]:Delta [from -0.5, 0.5]@FOOT%d E GM'ed, requested Q=%d", ifoot, Q_target),
		delta_binning, foot_binning
	};
	auto h2_e_vs_d = TH2P {
		Form("((h2_f))Cluster energy [ADC]:Delta [from -0.5, 0.5]@FOOT%d E fully calib, requested Q=%d", ifoot, Q_target),
		delta_binning, foot_binning
	};

	/* Main histogram, where we cut on δ */
	auto h2_e0_vs_x_mid = TH2P {
		Form("((h2_rc))Cluster raw energy [ADC]:Strip number [0..640]@FOOT%d central strip hit (dD=%.2f), requested Q=%d", ifoot, delta_cut, Q_target),
		bins_per_asic*TFOOTCalCont::N_ASIC, 0, RNFOOTCal::N_STRIPS, foot_binning
	};
	auto h1_e0_mid = TH1P {
		Form("((h1_e_mid))Cluster raw energy [ADC]@FOOT%d central strip hit (dD=%.2f), requested Q=%d", ifoot, delta_cut, Q_target),
		0xB2FD30_c, foot_binning
	};

    DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
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
                option::PostfixText{mnd::msg("[%u] (Re)gain FOOT match analysis (per event: %s)", i, fname.c_str())},
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
            if(mnd::PrintProgress(bar, entryId, nentries, 1000))
				progress.print_progress();

            const auto& sci21 = frs->cal.sci[0];
            const auto& sci22 = frs->cal.sci[1];
            const auto& sci31 = frs->cal.sci[2];

            if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
            if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
            if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);

            if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
            if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
            if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			bool has_valid_track = false;

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
				const double e_gm = e0 * gm;
				const double e = e_gm / foot_param->de.CorrectionFactor(delta);

				h1_delta.Fill(delta);

				h2_e0_vs_x.Fill(cx, e0);
				h2_egm_vs_d.Fill(delta, e_gm);
				h2_e_vs_d.Fill(delta, e);

                if(mnd::IsInside(delta, delta_interval)) {
					h2_e0_vs_x_mid.FillInside(cx, e0);
					h1_e0_mid.Fill(e0);
				}
				
				h1_foot_qc.Fill(t.Q);
				has_valid_track = true;
            } // for(const auto& t : foot->track)
	
			if(has_valid_track) {
				if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
				if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
				if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);
			}
        }
		if(!bar.is_completed())
			bar.mark_as_completed();
    });
    show_console_cursor(true);

	/* Idea is the following. Gain isn't always the same,.. some strips require higher gain for
	 * lower values. Simply to line up the total cluster energy values to the average ADC, across the detector. */
	
	/* Do the small fit in the 1D plot. */	
	TH1D* const h = &h1_e0_mid.h;
	auto [fitr_mid, err_mid] = GaussFitMax(h, sratio, niter);
	WARN("[CENTRAL HIT] 1D projection yields: max: %.2f, gauss fit max (around this max+-%.1f sigma): %.2f +- %.2f\n",
		h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ), sratio,
		fitr_mid[1], err_mid[1]);
	
	const double mean_mid = fitr_mid[1];

	std::vector<TGraph*>      profile_fit, gauss_fit;
	std::vector<TGraphErrors*> profile_raw, gauss_raw;

	const std::vector<int> fit_these_asics = {2,3,4,5,6,7};
	auto contains = [](const auto& v, const typename std::decay_t<decltype(v)>::value_type& val) -> bool {
		return std::find(v.begin(), v.end(), val) != v.end();
	};

	constexpr static int N_NEEDED_ENTRIES = 400;
	constexpr static int N_LOWEST_ENTRIES = 10;

	/* Try to fit a spline(s) for some given ASICs. */
	if(do_fit.is_some()) {
		FOOTGainParam pp = foot_param->gain;

		const auto& v = do_fit.unwrap();
		std::vector<double> fit_params;

		for(int a=0; a < TFOOTMapCont::N_ASIC; ++a) {
			FOOTAsicGainParam& asic = pp.fit.at(a);
			
			FMultiPoly* mp = asic.GetPoly(Q_target);
			if(!mp) {
				asic.multi_poly.emplace_back(Q_target);
				mp = &asic.multi_poly.back();
			}
			
			double x_lo  = (a) * 64 + 0.00001;
			double x_hi = (a+1) * 64 - 0.00001;
			
			if( contains(v, a) ) {
				auto [rg, graw, gfit] = FitSpline<fit_info::GAUSS_MAX> (
					poly_deg, h2_e0_vs_x_mid, x_lo, x_hi, 40, sratio, niter /*, Verbosity::CHATTY */
				);
				auto [rp, praw, pfit] = FitSpline<fit_info::PROFILE_MAX> (
					poly_deg, h2_e0_vs_x_mid, x_lo, x_hi, 40 /*, whatever, Verbosity::CHATTY */
				);
				gauss_fit.push_back(gfit);
				gauss_raw.push_back(graw);
				profile_fit.push_back(pfit);
				profile_raw.push_back(praw);
				
				if(take == Take::gauss or take == Take::gauss_fit_only) {
					mp->pol = std::vector<double>(rg.begin(), rg.end());
				} else {
					mp->pol = std::vector<double>(rp.begin(), rp.end());
				}
			}
			else {
				TAxis *xax = h2_e0_vs_x_mid->GetXaxis();
				int firstbin = xax->FindBin(x_lo);
				int lastbin = xax->FindBin(x_hi);

				auto pasic = std::unique_ptr<TH1D>(h2_e0_vs_x_mid->ProjectionY("__py", firstbin, lastbin));
				pasic->SetDirectory(nullptr);

				double profile_mean, gauss_mean;
				if(pasic->Integral() >= N_NEEDED_ENTRIES) { /* If it contains more than 500 events, we can sample it. */
					profile_mean = pasic->GetXaxis()->GetBinCenter( pasic->GetMaximumBin() );
					auto [pg0, err_pg0] = GaussFitMax( pasic.get(), sratio );
					
					if(!std::isfinite(pg0[1])) {
						WARN("ASIC[%d]: GaussFitMax produced non-finite mean value? Setting it to mean_mid.\n", a);
						gauss_mean = mean_mid;
					} else {
						gauss_mean = pg0[1];
					}
				} else { /* No clue. Just take profile mean the mean, but gauss is invalidated. */
					profile_mean = pasic->GetMean();
					gauss_mean = mean_mid;
				}
				TGraph* gfit = new TGraph(60);
				TGraph* pfit = new TGraph(60);
				for(int i=0; i<60; ++i) {
					double x = x_lo + (i+0.00001) * (x_hi - x_lo)/59;
					gfit->SetPoint(i, x, gauss_mean);
					pfit->SetPoint(i, x, profile_mean);
				}
				if(pasic->Integral() >= N_NEEDED_ENTRIES) {
					gfit->SetLineColor(gCol_); gfit->SetLineWidth(4);
					pfit->SetLineColor(pCol_); pfit->SetLineWidth(4);
				} else {
					gfit->SetLineColor(gCol_ + 1); gfit->SetLineWidth(12);
					pfit->SetLineColor(pCol_ + 1); pfit->SetLineWidth(12);
				}
				gauss_fit.push_back(gfit);
				profile_fit.push_back(pfit);

				mp->pol = std::vector<double>(1);
				if(pasic->Integral() < N_NEEDED_ENTRIES) {
					if(take == Take::gauss) { mp->pol[0] = gauss_mean; }
					else if(pasic->Integral() >= N_LOWEST_ENTRIES) { mp->pol[0] = profile_mean; }
					else { mp->pol[0] = gauss_mean; }; // in case if there are really no entries there...
				}
				else { // integral >= N_NEEDED_ENTRIES
					mp->pol[0] = gauss_mean;
				}
			}

			std::sort( asic.multi_poly.begin(), asic.multi_poly.end() );
		}
		auto j = nlohmann::json(pp);
		if(add_asic_indicator) {
			auto& fits = j.at("fit");
			int asic_i=0;
			for(auto& fit : fits) {
				fit["asic"] = asic_i++;
			}
		}

		std::cout << "\"gain\": " << j.dump(4) << std::endl;
		WARN("Average value: %.5f (bin-center) and %.5f (gauss-fit-center)\n",
			h1_e0_mid->GetBinCenter(h1_e0_mid->GetMaximumBin()), mean_mid);
	}

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, foot_binning[1],
				                i * 64, foot_binning[2]);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}
	
	[[maybe_unused]] TCanvas *c = new TCanvas("RawFOOT", Form("FOOT%d central", ifoot), 2400, 1400);
	h2_e0_vs_x_mid->SetStats(false);
	h2_e0_vs_x_mid->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(h2_e0_vs_x_mid->GetYaxis()->GetXmin());
		l->SetY2(h2_e0_vs_x_mid->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	for(auto* pfit : profile_fit) pfit->Draw("L SAME");
	for(auto* praw : profile_raw) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw) graw->Draw("P SAME");
	
	auto l = new TLegend(0.1,0.75,0.4,0.9);
	l->AddEntry(h2_e0_vs_x_mid, "Non-gain matched cluster energy (central hits)");
	if(!gauss_fit.empty())
		l->AddEntry(gauss_fit.front(), Form("Gaussian fit +-%.1f sigma around peak", sratio));
	if(!profile_fit.empty())
		l->AddEntry(profile_fit.front(), "TProfile fit");
	if(show_old.is_some()) {
		const auto& gain = foot_param->gain;
		const std::vector<i32>& Zs = show_old.unwrap();
		srand(time(NULL));
		const u32 v0 = (u32)rand();
		for(i32 Z : Zs) {
			try {
				auto [g,_] = gain.GetRefZGraph(Z);
				g->Draw("L SAME");
				g->SetLineColor( mnd::col::Col(v0 + Z) );
				g->SetLineStyle(kDashed);
				g->SetLineWidth(8);
				l->AddEntry(g, Form("Current gain curve (Z=%d) from the setup file", Z));
			} catch(std::exception const& e) {
				WARN("%s\n", e.what());
			}
		}
	}
	gStyle->SetLegendTextSize(0.021); l->Draw();
	gPad->GetCanvas()->ToggleEventStatus();

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21.Draw();
	cs->cd(2); h1_sci22.Draw();
	cs->cd(3); h1_sci31.Draw();
	cs->cd(4); h1_sci21_cut.Draw();
	cs->cd(5); h1_sci22_cut.Draw();
	cs->cd(6); h1_sci31_cut.Draw();

	TCanvas* cd = new TCanvas("delta_energy", "Delta & 1D energy", 1400, 800);
	cd->Divide(3,2);
	cd->cd(1); h1_delta.Draw();
	cd->cd(2); h2_egm_vs_d.Draw("COLZ"); gPad->SetLogz();
	cd->cd(5); h2_e_vs_d.Draw("COLZ"); gPad->SetLogz();

	cd->cd(4); h2_e0_vs_x.Draw("COLZ"); gPad->SetLogz();
	cd->cd(3); h1_e0_mid.DrawAndFit(sratio, kRed+1, 3/* line width */, niter);
	cd->cd(6);
	new PLatex(0.08,
		Form("Requested charge: Q=%d", Q_target),
		Form("Detector ID: FOOT%d", ifoot),
		Form("Delta cut: #pm%.3f", delta_cut),
		(take == Take::gauss) ? Form("Gauss fits size ratio: %.2f sigma", sratio)
			: "Fit taken from profile (violet curve)",
		Form("Mean mid: %.2f", (take == Take::gauss or take == Take::gauss_fit_only) ? mean_mid
			: h->GetXaxis()->GetBinCenter(h->GetMaximumBin()))
	);

	TCanvas *cf = new TCanvas("ft", "FOOT Track Selection", 1400,800);
	cf->Divide(2,1);
	cf->cd(1); h1_foot_q.Draw();
	cf->cd(2); h1_foot_qc.Draw();

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOT%d", ifoot),
		Form("Z_%d",  Q_target),
		mnd::rootfile_names_concatenated(fileNames)
	});

	WARN("End-of-main");
	rootApp.Run();
#else
    WARN("\n" BOLD ".. MND_FOOTTRACK_DEBUG not compiled in. This program does nothing.\n" KNRM);
#endif
    return 0;
}
