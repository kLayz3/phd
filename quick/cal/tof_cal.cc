/* Calibrate the time-of-flight between:
 * [0] : SCI21 -> SCI22 , to indicate "possible" velocity at S2
 * [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity at S3.
 *
 * No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run.
 * All of these measurements must be done with 12C files. There are 4 of relevance. */

#include "monad/monad.hxx"

#include "TCanvas.h"
#include "util/CLI.h"

#include "util/GaussFitMax.hxx"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/PolyFitter.h"
#include "util/Geometry.h"
#include "util/FitDrawer.hxx"
#include "util/Tracking.h"
#include "util/RunsheetParser.h"
#include "util/MPhysics.h"

#include "common/MacroCommon.hxx"

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
    CLI::App app{"Calibrate the time-of-flight between:\n\
                  [0] : SCI21 -> SCI22 , to indicate \"possible\" velocity at S2\n\
                  [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity and A/Q at S3.\n\
                  No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run (due to thick target).\n\
                  All of these measurements must be done with 12C files. There are 3 or 4 files of relevance."};

	std::vector<std::string> f;
	std::vector<TPCRef> ref{};
	std::array<double,3> dt_cut22_31 = {1000, -100, 100};
	std::array<double,3> dt_cut21_22 = {1000, -100, 100};
    double sratio = GAUSS_FIT_SIDE_RATIO_DEFAULT;
    double niter = 2;

	bool q = false;
	constexpr A3 qhist_df = {1000,0,4000};
	std::map<u32, A3> qhist_range {};
	unsigned short line_size = 4;
	std::vector<canvas::Extension> save = {};

    add_logged_option(app, "-f,--file", f, "Pass one or more file names and corresponding 2 brho's.")
		->delimiter(',');
    add_logged_option<DisplayDefault::No>(app, "-r, --ref", ref,
		"Select which TPC's (either with index: 0,1,2, or with a label: 21,22,23) make the reference. \
		Select by '0/1' which delay lines get included into the measurement. ")
		->type_name("[INT|LABEL:BOOL,BOOL;...]")
		->delimiter(';');
	add_logged_option(app, "--dt-cut-22-31", dt_cut22_31,
        "Delta T (SCI31 - SCI21) cut, in TDC units [25ps]")
		->delimiter(',');
	add_logged_option(app, "--dt-cut-21-22", dt_cut21_22,
        "Delta T (SCI22 - SCI21) cut, in TDC units [25ps]")
		->delimiter(',');
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_logged_flag(app, "-q,--charge", q, "Use this opportunity to also fit the charge. "
		"Must have the proper QDC pedestals already in the file.");
	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
	add_logged_option(app, "--qhist-range", qhist_range,
		"QDC value histogram range for specific scintillators")
		->default_str(mnd::streamable(qhist_df));
	add_logged_option(app, "-o,--save", save, "Save the resulting canvases as one or more extensions.")
		->delimiter(',');

    bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(f.size() < 3) {
		WARN("To continue, must supply at least 3 file names!\n"); return 0;
	}
    if(ref.size() < 2)
		ERROR("At least two valid referent TPC's must be given.\n");
    for(const auto& tpc : ref) {
		if(!tpc) { // operator bool()
			std::cerr << tpc << std::endl;
			ERROR("TPC invalid. Must be 0,1,2 and at least one dl flagged as valid.");
		}
	}
	mnd::fs::load_runsheet();
	TApplication rootApp("app", 0, 0);

    const u32 SCI21_I = RNFRSCal::SCI21_I;
    const u32 SCI22_I = RNFRSCal::SCI22_I;
    const u32 SCI31_I = RNFRSCal::SCI31_I;
    const u32 SCI_S2_I[2] = { SCI21_I, SCI22_I };
    const size_t nfiles = f.size();
   
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
    std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
    {
		std::unique_ptr<TFile> fhandle = std::make_unique<TFile>(f.front().c_str(), "READ");
        get_obj(fhandle, tpc_params, "FRS_tpc_parameters");
        get_obj(fhandle, sci_params, "FRS_sci_parameters");
    }
    constexpr auto N_TPC = TPCParam::N_S2_TPC;
	const Arr2<double, N_TPC, 2> zDL = TFRSCalCont::z_s2_tpc_delay_lines(tpc_params);

	std::array<SCIDEIntoQConverter const*, RNFRSCal::N_VALID_SCI> qdc_cvt{{nullptr}};
	if(q) {
		/* Try to find the pedestal values from the SCI's given. */
		for(u32 s = 0; s<RNFRSCal::N_VALID_SCI; ++s) {
			const SCIParam& scip = sci_params->at(s);
			const std::vector<SCIDEIntoQConverter>& vcvt = scip.de_to_q;
			if(vcvt.size() == 0) {
				WARN("Requested SCI%s QDC vs. β calibration, but no QDC pedestals in the file found. "
					"Is ok - will be skipping this scintillator.\n", label[s]);
				continue;
			}
			const SCIDEIntoQConverter* cvt =
				  (scip.SetConverter(f.front()) != 0)
				? scip.GetConverter()
				: &vcvt.front();

			if( ! mnd::isfinite(cvt->pedestal.left, cvt->pedestal.right) ) {
				WARN("Requested SCI%s QDC vs. β calibration, found a converted object, but QDC pedestals in the file not valid? "
					"Is ok - will be skipping this scintillator.\n", label[s]);
				continue;
			}
			qdc_cvt[s] = cvt;
		}
	}

    struct HistCont {
        TH1P *h_22_31, *h_21_22, *h_theta;
		std::array<TH1P*, RNFRSCal::N_VALID_SCI> h_q;
        TH2P *h_track_x;
    };
    std::vector<HistCont> hist;

	std::vector<double> x0, y0; // fitting containers ToF [0]
	std::vector<double> x1, y1; // fitting containers ToF [1]
    u32 i_clb_pnt = 1;
	
	std::array<mnd::col::RGBA, RNFRSCal::N_VALID_SCI> q_col;
	std::array<
		std::vector<double>,
		RNFRSCal::N_VALID_SCI
	> beta; /* beta[s] -> {beta_1, beta_2, beta_3, ...} => velocities for each of the calib files. */

	for(size_t i = 0; i<q_col.size(); ++i) {
		q_col[i] = mnd::col::next_col();
	}

	const auto run_info_example = mnd::QueryRunsheet(f.front());
	for(const auto& fileName : f) {
		const auto run_info = mnd::QueryRunsheet(fileName);
		WARN("[%s] runsheet row parsed as: %s\n", fileName.c_str(), mnd::streamable(run_info).c_str());
		
		if( run_info.primary.A != run_info.secondary.A
		 || run_info.primary.Z != run_info.secondary.Z) {
			ERROR("File: '%s', but primary beam: %s and secondary beam: %s do not match? This does not qualify as a calib run!\n",
				fileName.c_str(), run_info.primary.chem_to_string().c_str(), run_info.secondary.chem_to_string().c_str());
		}
		if( run_info.primary   != run_info_example.primary
		 || run_info.secondary != run_info_example.secondary) {
			ERROR("File/Initial file: '%s'/'%s' mismatch between either primary  %s/%s or secondary beams: %s/%s ?\n",
				fileName.c_str(), f.front().c_str(),
				run_info.primary.chem_to_string().c_str(), run_info_example.primary.chem_to_string().c_str(),
				run_info.secondary.chem_to_string().c_str(), run_info_example.secondary.chem_to_string().c_str());
		}

		const double A = run_info.secondary.A;
		const double Z = run_info.secondary.Z;

        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
        const double beta_ta  = phy::Beta(A, Z, phy::Brho_t{run_info.brho[mnd::TA_S1]});
        const double beta_inc = phy::Beta(A, Z, phy::Brho_t{run_info.brho[mnd::S1_S2]});
        const double beta_out = phy::Beta(A, Z, phy::Brho_t{run_info.brho[mnd::S2_S3]});
        const double beta_s4  = phy::Beta(A, Z, phy::Brho_t{run_info.brho[mnd::S3_S4]});
		WARN("[%s] calculated:\n"
			"β(TA-S1): " BOLD "%.4f" KNRM " , "
			"β(S1-S2): " BOLD "%.4f" KNRM " , "
			"β(S2-S3): " BOLD "%.4f" KNRM " , "
			"β(S3-S4): " BOLD "%.4f" KNRM " \n", fileName.c_str(), beta_ta, beta_inc, beta_out, beta_s4);
		
		beta[0].push_back(beta_inc); // SCI21
		beta[1].push_back(beta_out); // SCI22
		beta[2].push_back(beta_out); // SCI31
		beta[3].push_back(beta_s4);  // SCI41

        /* Containers for TPC extrapolation. */
        std::vector<double> xe, ye, ze;

		auto* h1_dt22_31 = new TH1P(Form("((h1_dt%u_22_31))Delta t [25 ps]@SCI31 - SCI22, TOF Point [%u]",
            i_clb_pnt, i_clb_pnt), kMagenta+i_clb_pnt, dt_cut22_31[0], dt_cut22_31[1], dt_cut22_31[2]);
		auto* h1_dt21_22 = new TH1P(Form("((h1_dt%u_21_22))Delta t [25 ps]@SCI22 - SCI21, TOF Point [%u]",
            i_clb_pnt, i_clb_pnt), kMagenta+i_clb_pnt, dt_cut21_22[0], dt_cut21_22[1], dt_cut21_22[2]);
        auto* h1_s2_angle = new TH1P(Form("((h1_s2_a%u))S2 polar angle [mrad]@TPC reference, point %u", i_clb_pnt, i_clb_pnt),
            kGreen -2, 500, 0, 20);
        auto* h2_track_x = new TH2P (
            Form("((h2_track_x%u))Track density (X) [mm]:Depth z [mm]@S2 area, TOF Point [%u]", i_clb_pnt, i_clb_pnt),
            800, 0, RNFRSCal::S2_LENGTH, 800, -60, 60);

		std::array<TH1P*, RNFRSCal::N_VALID_SCI> h1_q {};
		for(i32 i=0; i<RNFRSCal::N_VALID_SCI; ++i) {
			const auto it = qhist_range.find(i);
			A3 range = (it != qhist_range.end())? it->second: qhist_df;
			h1_q[i] = new TH1P(Form("((h1_q_%d))SCI%s QDC mean@TOF Point [%u]", i, label[i], i_clb_pnt),
				q_col[i], range[0], range[1], range[2]);
		}

        ProgressBar bar {
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{":)"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg("ToF Calibration (per event: %s)", fileName.c_str())},
                option::ForegroundColor{Color::magenta},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
        };
        const size_t nentries = ntuple->GetNEntries();

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
		    mnd::PrintProgress(bar, entryId, nentries, 500);

            const auto& sci21 = frs->sci[SCI21_I];
            const auto& sci22 = frs->sci[SCI22_I];
			const auto& sci31 = frs->sci[SCI31_I];

            /* All ToF measuring stations must be single hit. No exception. */
			if(sci22.hits.size() != 1 or sci31.hits.size() != 1 or sci21.hits.size() != 1)
                continue;

            /* We want to measure the S2 angle,.. important for SCI21 - SCI22 ToF ? */
            xe.clear(); ye.clear(); ze.clear();

            /* Find the reference containers. */
            for(const auto& id : ref) {
                u32 i = id.n;
                const auto& tpc = frs->tpc[i];
                for(u32 d : {0,1}) {
                    if(!id.use[d] or tpc.hits[d].size() != 1)
                        continue;
                    const RNTPCCal::Measurement& hit = tpc.hits[d].front();
                    const double x = hit.X();
                    const double y = hit.Y();
                    if(!std::isfinite(x) or !std::isfinite(y))
                        continue;
                    xe.push_back( x );
                    ye.push_back( y );
                    ze.push_back( zDL[i][d] );
                }
            }
            if(xe.size() < 3 or ye.size() < 3) continue;
            const auto fx = PolyFit<1>(ze, xe);
            const auto fy = PolyFit<1>(ze, ye);
            Line3D s2_track{ fx, fy };
            const double theta = s2_track.GetSpherical().theta * 1000; // mrad
            h1_s2_angle->Fill(theta);
			double dt22_31 = sci31.hits[0].t - sci22.hits[0].t;
			double dt21_22 = sci22.hits[0].t - sci21.hits[0].t;
			h1_dt22_31->Fill(dt22_31);
			h1_dt21_22->Fill(dt21_22);
            FillTrack(*h2_track_x, fx);

			for(u32 s=0; s<RNFRSCal::N_VALID_SCI; ++s) {
				/* If q==false, these are all silent nullptrs. */
				if( qdc_cvt[s] == nullptr) continue;

				double e = qdc_cvt[s]->E( frs->sci[s] );
				h1_q[s]->Fill(e);
			}

        } // for(auto entryId : *ntuple)

        hist.push_back({
            .h_22_31 = h1_dt22_31,
            .h_21_22 = h1_dt21_22,
            .h_theta = h1_s2_angle,
			.h_q = std::move(h1_q),
            .h_track_x = h2_track_x
		});
        auto [res22_31, _ ] = GaussFitMax(*h1_dt22_31, sratio, niter);
        auto [res21_22, __] = GaussFitMax(*h1_dt21_22, sratio, niter);
        x0.push_back( res22_31[1] ); // gauss mean value
        y0.push_back( 1.0 / beta_out );

        /* For S21-S22 β use their average. Normally, average path should also depend on theta,
         * but since it's mostly straight, and correction is O(theta^2), ignore it. */
        x1.push_back( res21_22[1] );
        y1.push_back( 2.0 / (beta_inc + beta_out) );
		
		if(!bar.is_completed()) bar.mark_as_completed();
        ++i_clb_pnt;
    } // for(const auto& fileName : f)
    if(hist.size() != nfiles) ERROR("Mismatch histogram and nfiles sizes.. ?");

    std::vector<double> fit_result_22_31, fit_result_21_22;
    auto [gerr0, g0] = FitAndDraw(1, x0, y0, {}, fit_result_22_31);
    auto [gerr1, g1] = FitAndDraw(1, x1, y1, {}, fit_result_21_22);
    gerr0->SetMarkerColor(kBlue - 1); gerr0->SetTitle("ToF Sci22 -> Sci31");
    gerr0->GetXaxis()->SetTitle("Mean ToF [25 ps], with offset");
    gerr0->GetYaxis()->SetTitle("1.0 / #beta");
    gerr1->SetMarkerColor(kRed - 1);  gerr1->SetTitle("ToF Sci21 -> Sci22");
    gerr1->GetXaxis()->SetTitle("Mean ToF [25 ps], with offset");
    gerr1->GetYaxis()->SetTitle("1.0 / #beta");
    g0->SetLineColor(kBlue);
    g1->SetLineColor(kRed);

    WARN("The formula is: " EMPH(1/β = a*ΔT + b\n));
    WARN("ToF S22 - S31: " EBOLD(b = %.6f; a = %.6f [1/(25ps)]\n), fit_result_22_31[0], fit_result_22_31[1]);
    WARN("ToF S21 - S22: " EBOLD(b = %.6f; a = %.6f [1/(25ps)]\n), fit_result_21_22[0], fit_result_21_22[1]);
	FRSToFSingle p22_31, p21_22;
	p22_31.combo = { SCI22_I, SCI31_I };
	p22_31.par[0] = fit_result_22_31[0]; p22_31.par[1] = fit_result_22_31[1];

	p21_22.combo = { SCI21_I, SCI22_I };
	p21_22.par[0] = fit_result_21_22[0]; p21_22.par[1] = fit_result_21_22[1];

	FRSToFParam tparam{};
	tparam.ToF.push_back(p22_31);
	tparam.ToF.push_back(p21_22);

	std::cout << nlohmann::json(tparam).dump(4) << std::endl;

    TCanvas *c = new TCanvas("Fit", "Fit", 1400, 700);
    c->Divide(2,1);
    c->cd(1); gerr0->Draw("AP"); gerr0->GetYaxis()->SetTitleOffset(1.2); g0->Draw("L SAME"); gPad->SetGrid();
    c->cd(2); gerr1->Draw("AP"); gerr0->GetYaxis()->SetTitleOffset(1.2); g1->Draw("L SAME"); gPad->SetGrid();

    TCanvas *c_raw = new TCanvas("RawToF", "RawToF", 2050, 1400);
    c_raw->Divide(4, nfiles);

    for(size_t i=0; i<nfiles; ++i) {
        auto [h_22_31, h_21_22, h_theta, h_q, h2_track_x] = hist[i];
        c_raw->cd(4*i + 1); h_22_31->DrawAndFit(sratio, kGreen, line_size, niter);
        c_raw->cd(4*i + 2); h_21_22->DrawAndFit(sratio, kGreen, line_size, niter);
        c_raw->cd(4*i + 3); h_theta->Draw();
        c_raw->cd(4*i + 4); h2_track_x->Draw("COLZ"); gPad->SetLogz();

		using mnd::hist::vline;
        double r = 0.72;
        TLine* line;
        for(int sci: {0,1}) {
            line = vline(h2_track_x, sci_params->at(SCI_S2_I[sci]).z0, r);
            line = vline(h2_track_x, sci_params->at(SCI_S2_I[sci]).z0, r);
            line->SetLineColor(kRed);
            line->SetLineStyle(2);
            line->SetLineWidth(4);
            line->Draw("SAME");
        }
    }

	/* ================= OPTIONAL Q vs. β calculation ==================== */
	if(q) {
		std::cerr << std::endl;
		TCanvas *c_q = new TCanvas("QDC", "QDC After Pedestal Subtraction", 2050, 1400);
		c_q->Divide(RNFRSCal::N_VALID_SCI , nfiles);
		std::array<std::vector<double>, RNFRSCal::N_VALID_SCI> q_fit_qdc {};
		
		for(size_t i=0; i<nfiles; ++i) {
			auto h_q = hist[i].h_q;
			if(h_q.size() != RNFRSCal::N_VALID_SCI)
				ERROR("Sizes mismatched! ...");
			for(size_t s=0; s < h_q.size(); ++s) {
				if(qdc_cvt[s] == nullptr)
					continue;

				c_q->cd(s+1 + i*RNFRSCal::N_VALID_SCI);
				auto [resq, _] = h_q[s]->DrawAndFit(sratio, kGreen, line_size, niter);
				q_fit_qdc[s].push_back( resq[1] );
			}
		}

		TCanvas *c_qf = new TCanvas("QDCFit", "QDC vs. β Fit", 2050, 1400);
		c_qf->Divide(3,2);

		std::vector<double> fit_q{};

		u32 i_q_cnv = 1;
		for(u32 s=0; s<RNFRSCal::N_VALID_SCI; ++s) {
			if(qdc_cvt[s] == nullptr)
				continue;
			if(q_fit_qdc[s].size() != beta[s].size())
				ERROR("SCI%s: sizes mismatch between β's array: %zu and number of fit pts %zu! Nfiles=%zu "
					"QDC converter 0x%p, all fine\n",
					label[s], beta[s].size(), q_fit_qdc[s].size(), f.size(), qdc_cvt[s]);
	
			std::vector<double> x = mnd::map(beta[s], [](double b) -> double { return 1.0 / (b*b); });
			std::vector<double> const& y = q_fit_qdc[s];

			WARN("SCI%s: fitsize: %zu, %zu. It is:\n%s\n%s\n", label[s],
				x.size(), y.size(), mnd::streamable(x).c_str(), mnd::streamable(y).c_str());
			auto [gerrQ, gQ] = FitAndDraw(1, x, y, {}, fit_q);
			gerrQ->SetMarkerColor(kRed - 1);
			gerrQ->SetTitle(Form("SCI%s: QDC vs. 1/#beta^{2} for calibration runs", label[s]));
			gerrQ->GetXaxis()->SetTitle("1/#beta^{2}");
			gerrQ->GetYaxis()->SetTitle("QDC Mean");
			gQ->SetLineColor(kRed);
			c_qf->cd(i_q_cnv++); gerrQ->Draw("AP"); gQ->Draw("L SAME"); gPad->SetGrid();
			WARN("Calculated dependency for SCI%s QDC vs. β: %s%s%s, with last QDC avg value: %s%.5f%s\n",
				label[s], KBH_MAG, mnd::streamable(fit_q).c_str(), KNRM, KBH_MAG, q_fit_qdc[s].back(), KNRM);

			SCIPrimary sp {};
			sp.beta = beta[s].back();
			sp.fit  = std::array{fit_q[0], fit_q[1]};
			
			nlohmann::json j{};
			j[ SCIParam::name_4 ] = sp;
			std::cout << BOLD "\"" << label[s] << "\": " << j.dump(4) << KNRM << ",\n";
		}
	}

    std::time_t now = std::time(nullptr);
    std::tm* tm = std::localtime(&now);
    char buffer[28];
    std::strftime(buffer, sizeof buffer,
        "%Y-%m-%d_%H-%M-%S", tm
    );
	canvas::save_all<canvas::Exe>(save, { buffer });
	WARN("End-of-main");
	rootApp.Run(); return 0;
}
