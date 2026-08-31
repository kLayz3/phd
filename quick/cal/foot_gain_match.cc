#include "util/CLI.h"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/FitSpline.h"

#include "TStyle.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

enum class Take { gauss, profile, gauss_fit_only };
enum class HitType { central, side };

using ShowOld = mnd::Option<std::vector<i32>>;
using DoFit   = mnd::Option<std::vector<i32>>;
namespace fs = std::filesystem;

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

int main(int argc, char* argv[]) {
	CLI::App app{"Perform gain matching of a specific FOOT detector."};
	
	std::string fileName = "";
	int ifoot = 0;
	int bins_per_asic = 64;
	DoFit do_fit {DoFit::No};
	double sratio = 0.9;
	u32 niter = 2;
	A3 foot_binning = {1000, 4, 4000};
	HitType hit_type = HitType::central;
	double delta_cut = 0.05;
	uint32_t mult_cut = 1; // any multiplicity below that is disallowed 
	int Q_target = 6;
	size_t poly_deg = 4;
	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	auto save = canvas::Extension::nil;
	Take take = Take::gauss_fit_only;
	ShowOld show_old { ShowOld::No };

	add_logged_option(app, "-f,--file", fileName, "Pass a file name.")
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-i,--foot-id", ifoot, 
		"Select which FOOT detector.")
		->check(CLI::Range(0, N_FOOT_DETECTORS-1));
	add_logged_option(app, "-q,--charge", Q_target,
		"Select which target charge we're roughly gating upon.")
		->check(CLI::Range(1,6));
	add_logged_option(app, "-p,--poly", poly_deg,
		"Select polynomial degree to fit the specified ASIC's.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-d,--delta", delta_cut,
		"Select the delta cut which will be applied to gate on very central hits.")
		->check(CLI::Range(0.01, 0.2));
	add_logged_option(app, "-c,--hit-type", hit_type,
		"Select to gate either on very central hit or very lateral (side) one.");
	add_logged_option(app, "-n,--foot-binning",foot_binning, "FOOT ADC binning")
		->delimiter(',');
	add_logged_option(app, "-b,--bins",bins_per_asic, "Binning Per ASIC")
		->check([](const std::string& match) -> std::string {
			int n = -1;
			if(!mnd::parse(match, n) or (TFOOTCalCont::N_STRIPS % n != 0)) 
				return mnd::msg("Must be integer parsable and evenly divisible by %d", TFOOTCalCont::N_STRIPS);
			return "";
		});
	add_logged_option(app, "-l,--fit", do_fit, "Select which ASICs to fit.")
		->delimiter(',')
		->type_name("INT[,INT...]")
		->check(CLI::RangeOrEmpty(0,10));
	add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber); 
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-m,--mult-cut", mult_cut, 
		"Multiplicity cut. Any clusters with multiplicity below selected one will be discarded.")
		->check(CLI::PositiveNumber);
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option(app, "-t,--take", take, 
		"Which type of projection fit to take. `gauss_fit_only` will only use gaussian spline for well sampled data.");
	add_logged_option(app, "-s,--show-old", show_old, 
		"Sequence of Z charges (, sep) to overlay their current gain match curve on the canvas.")
		->delimiter(',')
		->type_name("INT[,INT...]");
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");
	
	bool test = false;
	bool test_py = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");
	add_logged_flag(app, "--test-py", test_py, "Test the Python Matplotlib renderer.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileName.length() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

	FOOTParam *foot_param;
	{
		auto f = std::make_unique<TFile>(fileName.c_str(), "READ");
		get_obj(f, foot_param, Form("FOOT%d_setup", ifoot));
	}
	
	TApplication rootApp("app", 0, 0);
	
	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", 0x009B2F_c, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, 500, 300, 4000);
	auto* h1_delta     = new TH1P("((h1_d0))Delta [from -0.5, 0.5]", kGreen-2, 150, -0.499, 0.499);
	auto* h1_delta_cut_mid = new TH1P("((h1_d1_mid))Delta [from -0.5, 0.5]", kRed-7, 150, -0.499, 0.499);

	auto* h1_foot_e_mid = new TH1P("((h1_e_mid))FOOT E [ADC units]@Central strip value", 0xB2FD30_c, (int)(1.5*foot_binning[0]), foot_binning[1], foot_binning[2]);

	auto* hit_energy_mid = new TH2P(Form("((h2_mid))Cluster energy [ADC]:Strip number [0..640]@FOOT%d Raw, Requested Q=%d", ifoot, Q_target),
		bins_per_asic*10, 0,640, foot_binning[0], foot_binning[1], foot_binning[2]);
	
	const size_t nentries = ntuple->GetNEntries();
	A2 delta_interval_1 { -delta_cut, delta_cut };
	A2 delta_interval_2 { -delta_cut, delta_cut };
	
	if(hit_type == HitType::side) {
		WARN("Requesting side hit. Fetching the delta information part first...\n");
		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			const auto& sci21 = frs->sci[0];
			const auto& sci22 = frs->sci[1];
			const auto& sci31 = frs->sci[2];
			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;
			
			for(const auto& cl : foot->fCl) {
				if(cl.fCM < mult_cut) continue;
				double delta = cl.Delta();
				h1_delta->Fill(delta);
			}
		}
		TAxis* ax = (*h1_delta)->GetXaxis();
		ax->SetRangeUser(0.2, 0.45);
		const double dhi = ax->GetBinCenter( (*h1_delta)->GetMaximumBin() );
		ax->SetRange(0,0); // unzoom
						   //
		ax->SetRangeUser(-0.45, -0.2);
		const double dlo = ax->GetBinCenter( (*h1_delta)->GetMaximumBin() );
		ax->SetRange(0,0); // unzoom
		delta_interval_1 = { dhi-delta_cut, dhi+delta_cut };
		delta_interval_2 = { dlo-delta_cut, dlo+delta_cut };
		(*h1_delta)->Reset("ICESM");
	}

	ProgressBar bar {
		option::BarWidth{50},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("Analysis (per event: %s)", fileName.c_str())},
			option::ForegroundColor{Color::green},
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector{FontStyle::bold}}
	};
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		mnd::PrintProgress(bar, entryId, nentries, 1000, mnd::dancer0, 0.35);

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
			if(cl.fCM < mult_cut) continue;

			double delta = cl.Delta();
			h1_delta->Fill(delta);
			
			int i = static_cast<int>( cl.fCX );
			double e = cl.fCE;

			if(mnd::IsInside(delta, delta_interval_1) or mnd::IsInside(delta, delta_interval_2)) {
				h1_delta_cut_mid->Fill(delta);
				hit_energy_mid->FillInside(i, e);
				h1_foot_e_mid->FillInside(e);
			}
		}
	}
	bar.mark_as_completed();
	show_console_cursor(true);

	/* Idea is the following. Gain isn't always the same,.. some strips require higher gain for
	 * lower values. Simply to line up the total cluster energy values to the average ADC, across the detector. */
	
	/* Do the small fit in the 1D plot. */	
	TH1D* const h = *h1_foot_e_mid;
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
					poly_deg, *hit_energy_mid, x_lo, x_hi, 40, sratio, niter /*, Verbosity::CHATTY */
				);
				auto [rp, praw, pfit] = FitSpline<fit_info::PROFILE_MAX> (
					poly_deg, *hit_energy_mid, x_lo, x_hi, 40 /*, whatever, Verbosity::CHATTY */
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
				TAxis *xax = (*hit_energy_mid)->GetXaxis();
				int firstbin = xax->FindBin(x_lo);
				int lastbin = xax->FindBin(x_hi);

				auto pasic = std::unique_ptr<TH1D>((*hit_energy_mid)->ProjectionY("__py", firstbin, lastbin));
				pasic->SetDirectory(nullptr);

				double profile_mean, gauss_mean;
				if(pasic->Integral() >= N_NEEDED_ENTRIES) { /* If it contains more than 500 events, we can sample it. */
					profile_mean = pasic->GetXaxis()->GetBinCenter( pasic->GetMaximumBin() );
					auto [pg0, err_pg0] = GaussFitMax( pasic.get(), sratio );
					gauss_mean = pg0[1];
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

		std::cout << "\"gain\": " << nlohmann::json(pp).dump(4) << std::endl;
		WARN("Average value: %.5f (bin-center) and %.5f (gauss-fit-center)\n",
			(*h1_foot_e_mid)->GetBinCenter((*h1_foot_e_mid)->GetMaximumBin()), mean_mid);
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
	hit_energy_mid->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1((*hit_energy_mid)->GetYaxis()->GetXmin());
		l->SetY2((*hit_energy_mid)->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	for(auto* pfit : profile_fit) pfit->Draw("L SAME");
	for(auto* praw : profile_raw) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw) graw->Draw("P SAME");
	
	auto l = new TLegend(0.1,0.75,0.4,0.9);
	l->AddEntry(*hit_energy_mid, "Non-gain matched cluster energy (central hits)");
	l->AddEntry(gauss_fit[0], Form("Gaussian fit +-%.1f sigma around peak", sratio));
	l->AddEntry(profile_fit[0], "TProfile fit");
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

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	const char* ht_info = (hit_type == HitType::central)? "central": "side";
	TCanvas* cd = new TCanvas("delta_energy", "Delta & 1D energy", 1400, 800);
	cd->Divide(2,2);
	cd->cd(1); h1_delta->Draw();
	cd->cd(2); h1_delta_cut_mid->Draw();
	cd->cd(3); h1_foot_e_mid->DrawAndFit(sratio, kRed+1, 3/* line width */, niter);
	cd->cd(4);
	new PLatex(0.08,
		Form("Requested charge: Q=%d", Q_target),
		Form("Detector ID: FOOT%d", ifoot),
		Form("Delta cut: #pm%.3f", delta_cut),
		Form("Hit type: %s", ht_info),
		Form("Cluster size >= %d", mult_cut),
		(take == Take::gauss) ? Form("Gauss fits size ratio: %.2f sigma", sratio)
			: "Fit taken from profile (violet curve)",
		Form("Mean mid: %.2f", (take == Take::gauss or take == Take::gauss_fit_only) ? mean_mid
			: h->GetXaxis()->GetBinCenter(h->GetMaximumBin()))
	);

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOT%d", ifoot),
		Form("Z_%d",  Q_target),
		std::filesystem::path(fileName).stem().string(),
		ht_info
	});

	if(test_py) {
		mnd::plot::Figure {}
			.plot(
				*h1_delta, mnd::plot::HistStyle{}
					.stairs()
					.label(R"($\delta$ distribution)")
					.fill()
					.line_width(2.2)
			).xlabel(R"($\delta\,[-0.5,\,0.5]$)")
			.ylabel(*h1_delta)
			.grid()
			.title(*h1_delta)
			.legend()
			.save(
				fs::path{"autosave"} / mnd::fs::current_executable_name() / "pysave.png"
			);
		mnd::plot::Figure {}
			.plot(
				*h1_delta, mnd::plot::HistStyle{}
					.stairs()
					.label(R"($\delta$ distribution)")
					.line_width(2.2)
			).xlabel(R"($\delta\,[-0.5,\,0.5]$)")
			.ylabel(*h1_delta)
			.grid()
			.title("Delta, non-fill 'gram")
			.legend()
			.save(
				fs::path{"autosave"} / mnd::fs::current_executable_name() / "pysave_nf.png"
			);
		mnd::plot::Figure {}
			.plot(
				*h1_delta, mnd::plot::HistStyle{}
					.points()
					.label(R"($\delta$ distribution)")
					.line_width(2.2)
			).xlabel(R"($\delta\,[-0.5,\,0.5]$)")
			.ylabel("Count")
			.grid()
			.title("Delta, points 'gram")
			.legend()
			.save(
				fs::path{"autosave"} / mnd::fs::current_executable_name() / "pysave_pt.png"
			);
	}

	WARN("End-of-main\n");

	rootApp.Run(); return 0;
}
