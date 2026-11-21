/* Calibrate the time-of-flight between:
 * [0] : SCI21 -> SCI22 , to indicate "possible" velocity at S2
 * [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity at S3.
 *
 * No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run.
 * All of these measurements must be done with 12C files. There are 4 of relevance. */

#include "cli/CLI11.hpp"
#include "monad/monad.hxx"

#include "TCanvas.h"
#include "util/CLI.h"

#include "util/GaussFitMax.hxx"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/FitDrawer.hxx"
#include "util/RunsheetParser.h"
#include "util/MPhysics.h"

#include "TApplication.h"

#include "TFRSCalCont.h"
#include "TFRSHitCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::geom;
using namespace mnd::col::literals;

const auto& label = RNFRSCal::sci_label;

int main(int argc, char* argv[]) {
	static_assert(RNFRSCal::N_VALID_SCI > 1,
		"Must have at least two or more valid scintillators compiled in!");

    CLI::App app{
		"Calibrate the time-of-flight between a pair of two scintillators selected by the -s,--sci flags.\n\
		All of these measurements must be done with 12C files. There are 3 or 4 files of relevance.\n\
		This program can also used to do β-calibration of the mean QDC values, based on the measured \n\
		velocity for both of the scintillators."
	};

	std::vector<std::string> f;
	A3 dt_cut = {5000, -1000, 1000};
	std::array<int,2> sci_i;
    double sratio = GAUSS_FIT_SIDE_RATIO_DEFAULT;
    double niter = 2;

	bool q = false;
	constexpr A3 qhist_df = {1000,0,4000};
	std::map<u32, A3> qhist_range {};
	unsigned short line_size = 4;
	std::vector<canvas::Extension> save = {};

	std::vector<u32> t_exclude;
	std::map<u32, std::vector<u32>> q_exclude;

    add_logged_option(app, "-f,--file", f, "Pass one or more file names and corresponding 2 brho's.")
		->delimiter(',');
	add_logged_option(app, "-s,--sci", sci_i, "Pass two scintillator indices.")
		->delimiter(',')
		->check(CLI::Range(0, RNFRSCal::N_VALID_SCI-1))
		->required();
    add_logged_option(app, "--dt-cut", dt_cut,
        "∆T (SCI[1] - SCI[0]) cut, in TDC units [25ps]")
		->delimiter(',');
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_logged_flag(app, "-q,--charge,!--no-charge", q, "Use this opportunity to also fit the average QDC value vs. β. "
		"Must have the proper QDC pedestals already in the file. Check cal/sci_qdc_ped program for that.");
	add_logged_option(app, "--q-exclude", q_exclude, "Exclude some measurement points (β¡, E₀¡), i∈{0,1,2...} "
		"(0-indexed) from the charge-analysis. Map entry should be either [0] or [1] denoting the first "
		"or the second scintillator, respectively. Map value is a sequence of points, ',' separated. "
		"If requested outside, is ignored.");
	add_logged_option(app, "--t-exclude", t_exclude, "Exclude the measurement point (β¡, ∆T¡), i∈{0,1,2...} "
		"(0-indexed) from the ToF analysis. If requested outside, is ignored.")
		->delimiter(',');

	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
	add_logged_option(app, "--qhist-range", qhist_range,
		"QDC value histogram range for all the scintillators. "
		"The two indices passed in -s,--sci should match some of the two inputs of the map. "
		"If not found, will assign the default value.")
		->default_str(fmt::format("{}", qhist_df));

	add_logged_option(app, "-o,--save", save, "Save the resulting canvases as one or more extensions.")
		->delimiter(',');

    bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(f.size() < 3) {
		ERROR("To continue, must supply at least 3 file names!\n");
	}
	
	if(sci_i[0] >= sci_i[1])
		ERROR("Must have SCI indices START < STOP.\n");

	mnd::fs::load_runsheet();
	TApplication rootApp("app", 0, 0);

    const size_t nfiles = f.size();

    std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
    {
		std::unique_ptr<TFile> fhandle = std::make_unique<TFile>(f.front().c_str(), "READ");
        get_obj(fhandle, sci_params, "FRS_sci_parameters");
    }

	std::array<SCIDEIntoQConverter,2> qdc_cvt {};
	if(q) {
		/* Try to find the pedestal values from the SCI's given. */
		for(int i : {0,1}) {
			const int s = sci_i[i];
			const SCIParam& scip = sci_params->at(s);
			const SCIDEIntoQConverter& cvt = scip.qdc;
			if(cvt.IsDefaulted() )
				ERROR("Requested SCI{} QDC vs. β calibration, but no valid QDC pedestals in the file found? "
					"Find the pedestals first, using the `cal/sci_qdc_ped` program.\n", label.at(s));
		
			qdc_cvt[i] = cvt;
		}
	}

    struct HistCont {
        TH1P *h1_dt;
		std::array<TH1P*, 2> h1_q;
		std::array<TH2P*, 2> h2_q_lr;
    };
    std::vector<HistCont> hist;
	hist.reserve(nfiles);

	std::vector<double> x, y; // fitting containers for Time-of-Flight
    u32 i_clb_pnt = 1;
	
	constexpr static std::array<mnd::col::RGBA, 2> q_col
		= { 0xAB0732_c, 0x01B7DB_c};

	std::vector<double> beta; // {beta_1, beta_2, beta_3, ...} => velocities for each of the calib files.
	beta.reserve(nfiles);

	const auto run_info_example = mnd::QueryRunsheet(f.front());
	for(const auto& fileName : f) {
		const auto run_info = mnd::QueryRunsheet(fileName);
		WARN("[{}] runsheet row parsed as: {}\n", fileName, run_info);
		
		if( run_info.primary.A != run_info.secondary.A
		 || run_info.primary.Z != run_info.secondary.Z) {
			ERROR("File: '{}', but primary beam: {} and secondary beam: {} do not match? This does not qualify as a calib run!\n",
				fileName, run_info.primary.chem_to_string(), run_info.secondary.chem_to_string());
		}
		if( run_info.primary   != run_info_example.primary
		 || run_info.secondary != run_info_example.secondary) {
			ERROR("File/Initial file: '{}'/'{}' mismatch between either primary  {}/{} or secondary beams: {}/{} ?\n",
				fileName, f.front().c_str(),
				run_info.primary.chem_to_string(), run_info_example.primary.chem_to_string(),
				run_info.secondary.chem_to_string(), run_info_example.secondary.chem_to_string());
		}

		const phy::Nucleus& ion = run_info.secondary;
		WARN("Secondary ion parsed as: {}\n", styled(ion.chem_to_string(), mnd::style::bold_magenta));

        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
        const double beta_ta  = phy::Beta(ion, phy::Brho_t{run_info.brho[mnd::TA_S1]});
        const double beta_s2_inc = phy::Beta(ion, phy::Brho_t{run_info.brho[mnd::S1_S2]});
        const double beta_s2_out = phy::Beta(ion, phy::Brho_t{run_info.brho[mnd::S2_S3]});
        const double beta_s4  = phy::Beta(ion, phy::Brho_t{run_info.brho[mnd::S3_S4]});
		WARN("[{}] calculated:\n"
			"β(TA-S1): {:.4f}, β(S1-S2): {:.4f}, β(S2-S3): {:.4f}, β(S3-S4): {:.4f}\n",
			fileName, styled(beta_ta, fmt::emphasis::bold), styled(beta_s2_inc, fmt::emphasis::bold),
			styled(beta_s2_out, fmt::emphasis::bold), styled(beta_s4, fmt::emphasis::bold));

		switch(sci_i[1]) {
			case 0: beta.push_back(beta_s2_inc); break; // SCI21
			case 1: beta.push_back(beta_s2_out); break; // SCI22
			case 2: beta.push_back(beta_s2_out); break; // SCI31
			case 3: beta.push_back(beta_s4);     break; // SCI41
			default:
				ERROR("Case not handled? Is: {}\n", sci_i[1]);
		}
		WARN("SCI{}-{} nominal β-identified as: {:.4f}\n",
			label.at(sci_i[0]), label.at(sci_i[1]), styled(beta.back(), mnd::style::bold_magenta));

		auto* h1_dt = new TH1P(Form("((h1_dt_%d))#Delta t [25 ps]@SCI%s - SCI%s, TOF Point [%u]",
            i_clb_pnt, label.at(sci_i[0]), label.at(sci_i[1]), i_clb_pnt),
			kMagenta+i_clb_pnt, dt_cut[0], dt_cut[1], dt_cut[2]);

		std::array<TH1P*, 2> h1_q {};
		std::array<TH2P*, 2> h2_q_lr {};
		for(int i : {0,1}) {
			const int s = sci_i[i];
			const auto it = qhist_range.find(s);
			A3 range = (it != qhist_range.end())? it->second: qhist_df;
			
			h1_q[i] = new TH1P(Form("((h1_q_%d))SCI%s QDC mean [QDC units]@TOF Point [%u]", i, label.at(s), i_clb_pnt),
				q_col[i], range[0], range[1], range[2]);

			h2_q_lr[i] = new TH2P(Form("((h2_q_%d))QDC left [QDC units]:QDC right [QDC units]@SCI%s, TOF Point [%u]", i, label.at(s), i_clb_pnt),
				range[0], range[1], range[2], range[0], range[1], range[2]);
		}

        const size_t nentries = ntuple->GetNEntries();
        ProgressBar bar {
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{":)"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{fmt::format("ToF Calibration {} (per event: {})", nentries, fileName)},
                option::ForegroundColor{Color::magenta},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
        };
	
		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
		    mnd::PrintProgress(bar, entryId, nentries, 500);

            const auto& sci0 = frs->sci.at( sci_i[0] );
            const auto& sci1 = frs->sci.at( sci_i[1] );

            /* All ToF measuring stations must be single hit. No exception. */
			if(sci1.hits.size() != 1 || sci0.hits.size() != 1)
                continue;

            const double dt = sci1.hits[0].t - sci0.hits[0].t;
			h1_dt->Fill(dt);

			if(q) {
				double e_sci0_l = qdc_cvt[0].E0_l( sci0 );
				double e_sci0_r = qdc_cvt[0].E0_r( sci0 );
				double e_sci1_l = qdc_cvt[1].E0_l( sci1 );
				double e_sci1_r = qdc_cvt[1].E0_r( sci1 );

				double e_sci0 = qdc_cvt[0].E0( sci0 );
				double e_sci1 = qdc_cvt[1].E0( sci1 );
				
				/* Problem with RAW QDC values is that it usually has binning issues,
				 * so add a small width to it. */
				h1_q[0]->Fill(e_sci0);
				h1_q[1]->Fill(e_sci1);
				h2_q_lr[0]->Fill(e_sci0_r, e_sci0_l);
				h2_q_lr[1]->Fill(e_sci1_r, e_sci1_l);
			}

        } // for(auto entryId : *ntuple)
	
        auto [result, _ ] = GaussFitMax(*h1_dt, sratio, niter);
        x.push_back( result[1] ); // gauss mean value
        y.push_back( 1.0 / beta.back() );

        hist.push_back({
            .h1_dt = h1_dt,
			.h1_q = std::move(h1_q),
			.h2_q_lr = std::move(h2_q_lr),
		});

		if(!bar.is_completed())
			bar.mark_as_completed();
		++i_clb_pnt;
    } // for(const auto& fileName : f)
	
    if(hist.size() != nfiles)
		ERROR("Mismatch histogram and nfiles sizes.. ?");

    std::vector<double> fit_result;
    auto [gerr, g] = FitAndDraw(1, x, y, {}, t_exclude, fit_result);
    gerr->SetMarkerColor(kBlue - 1);
	gerr->SetTitle(Form("ToF SCI%s -> SCI%s", label.at(sci_i[0]), label.at(sci_i[1])));
    gerr->GetXaxis()->SetTitle("Mean ToF [25 ps], with offset");
    gerr->GetYaxis()->SetTitle("1.0 / #beta");
    g->SetLineColor(kBlue);

	WARN("The formula is: {}\n", styled("1/β = a*ΔT + b", mnd::style::bold_yellow));
	WARN("Time-of-Flight SCI{} - SCI{}: {}\n", label.at(sci_i[0]), label.at(sci_i[1]),
		styled(fmt::format("b = {:.6f}; a = {:.6f} [1/(25ps)]", fit_result[0], fit_result[1]),
			fmt::emphasis::bold));
	FRSToFSingle p{};
	p.combo = { (u32)sci_i[0], (u32)sci_i[1] };
	p.par[0] = fit_result[0];
	p.par[1] = fit_result[1];

	fmt::println("{}", nlohmann::json(p).dump(4));

    TCanvas *c = new TCanvas("ToF", "Time-of-Flight to do β-calculation", 2150, 1250);
    c->Divide(2,1);
    c->cd(1); gerr->Draw("AP"); gerr->GetYaxis()->SetTitleOffset(1.2); g->Draw("L SAME"); gPad->SetGrid();

    TCanvas *c_raw = new TCanvas("RawToF", "RawToF", 2050, 1400);
    c_raw->Divide(2, nfiles);

    for(size_t i=0; i<nfiles; ++i) {
        auto [h1_dt, _, __] = hist[i];
        c_raw->cd(2*i + 1); h1_dt->DrawAndFit(sratio, kGreen, line_size, niter);
	}

	/* ================= OPTIONAL Q vs. β calculation ==================== */
	if(q) {
		fputc('\n', stderr);
		TCanvas *c_q = new TCanvas("QDC", "QDC After Pedestal Subtraction", 2050, 1400);
		c_q->Divide(4 , nfiles);

		std::array<
			std::vector<double>, 2
		> q_fit_qdc {};
		
		for(size_t n=0; n<nfiles; ++n) {
			auto h_q    = hist[n].h1_q;
			auto h_q_lr = hist[n].h2_q_lr;
			if(h_q.size() != 2)
				ERROR("Sizes mismatched! ...");
			for(auto i: {0,1}) {
				c_q->cd(2*i+1 + 4*n);
				auto [result_q, _] = h_q.at(i)->DrawAndFit(sratio, kGreen, line_size, niter);
				q_fit_qdc[i].push_back( result_q[1] );

				c_q->cd(2*i+2 + 4*n);
				h_q_lr[i]->Draw("COLZ");
			}
		}

		TCanvas *c_qf = new TCanvas("QDCFit", "QDC vs. 1/β² fit", 2050, 1400);
		c_qf->Divide(2,1);

		std::vector<double> fit_q{};

		for(int i : {0,1}) {
			const int s = sci_i[i];

			if(q_fit_qdc[i].size() != beta.size())
				ERROR("SCI{}: sizes mismatch between β's array: {} and number of fit pts {}! Nfiles={} "
					"QDC converter {}, all fine\n",
					label[s], beta.size(), q_fit_qdc[i].size(), f.size(), qdc_cvt[i] );
	
			std::vector<double> x = mnd::map(beta, [](double b) -> double { return 1.0 / (b*b); });
			std::vector<double> const& y = q_fit_qdc[i];

			WARN("SCI{}: fitsize: {}, {}. It is:\n{}\n{}\n", label[s],
				x.size(), y.size(), x, y);
			auto [gerrQ, gQ] = FitAndDraw(1, x, y, {}, q_exclude[i], fit_q);
			gerrQ->SetMarkerColor(kRed - 1);
			gerrQ->SetTitle(Form("SCI%s: QDC vs. 1/#beta^{2} for calibration runs", label[s]));
			gerrQ->GetXaxis()->SetTitle("1/#beta^{2}");
			gerrQ->GetYaxis()->SetTitle("QDC Mean");
			gQ->SetLineColor(kRed);
			c_qf->cd(i+1); gerrQ->Draw("AP"); gQ->Draw("L SAME"); gPad->SetGrid();
			WARN("Calculated dependency for SCI{} QDC vs. β: {}, with last QDC avg value: {:.5f}\n",
				label[s], styled(fit_q, mnd::style::bold_magenta), styled(q_fit_qdc[i].back(), mnd::style::bold_magenta));

			SCIPrimary sp {};
			sp.beta = beta.back();
			sp.fit  = std::array{fit_q[0], fit_q[1]};
			
			nlohmann::json j{};
			j[ SCIParam::name_4 ] = sp;
			fmt::println("{},", styled(fmt::format("\"{}\": {}", label[s], j.dump(4)), fmt::emphasis::bold));
		}
	}

	const std::string info = std::string{"tof"}
		+ fmt::format("{}-{}{}", label.at(sci_i[0]), label.at(sci_i[1]), q ? "-q" : "");
	canvas::save_all<canvas::Exe>(save, { info });
	WARN("End-of-main");
	rootApp.Run(); return 0;
}
