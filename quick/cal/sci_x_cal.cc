#include "util/CLI.h"

#include "util/GaussFitMax.hxx"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.hxx"
#include "util/PolyFitter.h"
#include "util/Tracking.h"
#include "util/FitSpline.h"

#include "common/MacroCommon.hxx"

#include "TApplication.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

 int main(int argc, char* argv[]) {
	CLI::App app{"Calibrate the X measuring dimension (roughly) of either SCI21 and SCI22"};
    
    std::vector<std::string> fileName{};
	std::vector<TPCRef> ref{}; 
	u32 i_sci = 0;
    bool do_diff = false;
	A3 binning_ref = {100, -30, 30};
	A3 binning_sci = {100, -30, 30};
	A2 fit_range = {-5, 5};
    u32 niter = 2;
    double sratio = 1.4;
	auto save = canvas::Extension::nil;

    add_logged_option(app, "-f,--file", fileName, "Pass one or more file names from the ToF cal run.")
		->delimiter(',')
        ->check(CLI::ReadPermissions);
    add_logged_option(app, "-i,--sci", i_sci, "Scintillator index; 0 => SCI21, 1 => SCI22.")
        ->check(CLI::Range(0,1));
    add_logged_option<DisplayDefault::No>(app, "-r, --ref", ref, 
		"Select which TPC's (either with index: 0,1,2, or with a label: 21,22,23) make the reference. \
		Select by '0/1' which delay lines get included into the measurement. ")
		->type_name("[INT|LABEL:BOOL,BOOL;...]")
		->delimiter(';');

	add_logged_option(app, "-x,--binning-ref",binning_ref, "Binning X.")
		->delimiter(',');
	add_logged_option(app, "-s,--binning-sci",binning_sci, "Binning Y. If difference toggle given, then Y becomes the difference axis (delta axis).")
		->delimiter(',');
	add_logged_option(app, "-r,--fit-range",fit_range, "Fit range [mm], along reference X-axis, given by binning-ref option.")
		->delimiter(',');
    add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
        ->check(CLI::PositiveNumber);
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
        ->check(CLI::PositiveNumber);
	add_logged_flag(app, "-d,--diff", do_diff, "For the correlation plot, have y-axis be the difference (ref - sci), instead of the actual SCI measurement.");
    add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileName.size() == 0) {
		WARN("To continue, must supply at least one file name!\n"); return 0;
	}
    const size_t nfiles = fileName.size();

	if(ref.size() < 2) 
		ERROR("At least two valid referent TPC's must be given.\n");
    for(const auto& tpc : ref) {
		if(!tpc) { // operator bool() 
			std::cerr << tpc << std::endl; 
			ERROR("TPC invalid. Must be 0,1,2 and at least one dl flagged as valid."); 
		}
	}
    const auto& label = TFRSCalCont::sci_label;
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;

	TApplication rootApp("app", 0, 0);
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.front().c_str(), "READ");
        get_obj(f, tpc_params, "FRS_tpc_parameters");
        get_obj(f, sci_params, "FRS_sci_parameters");
    }

	const auto& sci_param = sci_params->at(i_sci);
	const double z0 = sci_param.z0;
	constexpr double TPC_WIDTH = 70.0;
	const Arr2<double, N_UPSTREAM_TPC, 2> zDL = [tpc_params](){
		Arr2<double,N_UPSTREAM_TPC, 2> z{};
		for(u32 i=0; i<N_UPSTREAM_TPC; ++i) {
			z[i][0] = tpc_params->at(i).z0 - TPC_WIDTH/4;
			z[i][1] = tpc_params->at(i).z0 + TPC_WIDTH/4;
		}
		return z;
	}();

    /* Containers for TPC extrapolation. */
	std::vector<double> xe, ye, ze;
    std::vector<TH1P*> g_h1_sci;
    std::vector<TH2P*> g_histd, g_h2_track_x, g_h2_track_y, g_h2_ab, g_h2_xy;
    u32 i_clb_pnt = 1;
    for(const auto& fname : fileName) {
        auto* h1_sci = new TH1P (
            Form("((h1_sci%u))SCI%s QDC mean [QDC units]@Calibration point %u", i_clb_pnt, label[i_sci], i_clb_pnt), 
            ORGB{0xCB00CB}, 500, 300, 4000
        );
        auto* histd = new TH2P(Form("SCI%s X -  TPC extr.[mm]:TPC extr. [mm]", label[i_sci]), 
            binning_ref[0], binning_ref[1], binning_ref[2],
            binning_sci[0], binning_sci[1], binning_sci[2]);

        auto* h2_track_x = new TH2P(Form("((h2_track_x%u))Track density (X) [mm]:Depth z [mm]@S2 area, clb point %u", i_clb_pnt, i_clb_pnt), 600, 0, z0 + 20, 600, -60, 60);
        auto* h2_track_y = new TH2P(Form("((h2_track_x%u))Track density (Y) [mm]:Depth z [mm]@S2 area, clb point %u", i_clb_pnt, i_clb_pnt), 600, 0, z0 + 20, 600, -60, 60);
        auto* h2_ab = new TH2P(Form("((h2_ab%u))Y-angle [mrad]:X-angle [mrad]@At SCI%s, clb point %u", i_clb_pnt, TFRSCalCont::sci_label[i_sci], i_clb_pnt),  100, -20, 20, 100, -20, 20);
        auto* h2_xy = new TH2P(Form("((h2_xy%u))Referent y-position [mm]:Referent X-position [mm]@At SCI%s, clb point %u", i_clb_pnt, label[i_sci], i_clb_pnt),
		    binning_ref[0],binning_ref[1],binning_ref[2],binning_ref[0],binning_ref[1],binning_ref[2]);

        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fname);
        ProgressBar bar {
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg("Analysis (per event: %s)", fname.c_str())},
                option::ForegroundColor{Color::red},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
        };
        const size_t nentries = ntuple->GetNEntries();

        for(auto entryId : *ntuple) {
            ntuple->LoadEntry(entryId);
		    mnd::PrintProgress(bar, entryId, nentries, 1000, mnd::as_span(mnd::dancer0), 0.35);

		    const auto& sci = frs->sci[i_sci];
            
            /* Care only about single-hit events. */
            if(sci.hits.size() != 1) continue;
            h1_sci->Fill(sci.E);
            const RNSciCal::Measurement& sci_hit = sci.hits.front();

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
                    xe.push_back(x);
                    ye.push_back(y);
                    ze.push_back(zDL[i][d]);
                }
            }
            if(xe.size() < 3 or ye.size() < 3) continue;
            const auto fx = PolyFit<1>(ze, xe);
            const auto fy = PolyFit<1>(ze, ye);
             
            /* Extrapolated positions at the FOOT: */
            const double xRef = fx[1] * z0 + fx[0];
            const double yRef = fy[1] * z0 + fy[0];

            /* SCI's just measure the x' ... as such should be invariant relative to y. 
             * but keep it anyway.. its w/e */
            FillTrack(*h2_track_x, fx);
            FillTrack(*h2_track_y, fy);
            h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);
            h2_xy->Fill(xRef, yRef);
            
            const double measurement = (do_diff)? (sci_hit.x - xRef): xRef;
            histd->Fill(xRef, measurement); 
        }
	    bar.mark_as_completed();
       
        constexpr int N_PTS_FOR_GRAPH = 30;
        auto [rg, gerr, g] = FitSpline<1, fit_info::GAUSS_MAX> (
            *histd, fit_range[0], fit_range[1], N_PTS_FOR_GRAPH, sratio, niter /*, Verbosity::CHATTY */
        );

        auto [offset, slope]  = rg;
        const double offset0 = sci_param.x_offset; // already inputted.
        const double slope0 = sci_param.x_factor;   // already inputted.

        char text0[20] = {'\0'};
        sprintf(text0, "[%u] For SCI%s:", i_clb_pnt, label[i_sci]);

        char text1[1024] = {'\0'};
        sprintf(text1, "[%u] Calculated: graph offset = %.5f, graph slope = %.5f.", i_clb_pnt, offset, slope);
        WARN("%s %s\n", text0, text1);

        char text2[1024] = {'\0'};
        sprintf(text2, "[%u] Currently in the setup file: (offset, slope) = (%.5f, %.5f)", i_clb_pnt, sci_param.x_offset, sci_param.x_factor);
        WARN("%s\n", text2);

        const double offset1 = (sci_param.x_offset - offset) / slope; // values to be written into file.
        const double slope1 = sci_param.x_factor / slope;             // values to be written into file.
        char result_slope[1024] = {'\0'};
        char result_offset[1024] = {'\0'};
        sprintf(result_offset, "x_offset: " KCYN "%.5f [mm]" KNRM, offset1);
        sprintf(result_slope, "x_factor: " KCYN "%.5f [mm/tdc]" KNRM, slope1);
        WARN("[%u] Recommended: %s\n%s\n", i_clb_pnt, result_offset, result_slope);
        WARN("[%u] ~~~~~ This is a relative change of: %.3f%% and %.3f%%\n", i_clb_pnt,
            100*std::abs((slope1-slope0)/slope0), 100*std::abs((offset1 - offset0)/offset0));
        
        g_histd.push_back(histd);
        g_h2_track_x.push_back(h2_track_x);
        g_h2_track_y.push_back(h2_track_y);
        g_h2_ab.push_back(h2_ab);
        g_h2_xy.push_back(h2_xy);
        g_h1_sci.push_back(h1_sci);
    
        ++i_clb_pnt;
    }
    assert(g_histd.size() == nfiles      && "Sanity #0");
    assert(g_h2_track_x.size() == nfiles && "Sanity #1");
    assert(g_h2_track_y.size() == nfiles && "Sanity #2");
    assert(g_h2_ab.size() == nfiles      && "Sanity #3");
    assert(g_h2_xy.size() == nfiles      && "Sanity #4");
    assert(g_h1_sci.size() == nfiles     && "Sanity #5");

	TCanvas *c = new TCanvas("TPC-Ref", "TPC-Ref", 2400, 1400);
    c->Divide(4, nfiles);
    for(size_t i=0; i<nfiles; ++i) {
        c->cd(4*i + 1); gPad->SetLogz(); g_h2_track_x[i]->Draw("COLZ");
        c->cd(4*i + 2); gPad->SetLogz(); g_h2_track_y[i]->Draw("COLZ");
        c->cd(4*i + 3); gPad->SetLogz(); g_h2_xy[i]->Draw("COLZ");
        c->cd(4*i + 3); gPad->SetLogz(); g_h2_ab[i]->Draw("COLZ");
    }
    
	TCanvas *cs = new TCanvas("SCI-E-Corr", "SCI-E-Corr", 2400, 1400);
    cs->Divide(2, nfiles);
    for(size_t i=0; i<nfiles; ++i) {
        cs->cd(2*i + 1); g_histd[i]->Draw("COLZ");
        cs->cd(2*i + 2); g_h1_sci[i]->Draw();
    }

    std::time_t now = std::time(nullptr);
    std::tm* tm = std::localtime(&now);
    char buffer[28];
    std::strftime(buffer, sizeof buffer,
        "%Y-%m-%d_%H-%M-%S", tm
    );
	canvas::save_all<canvas::Exe>(save, { Form("sci%s", label[i_sci]), buffer });

	WARN("End-of-main");
	rootApp.Run(); return 0;
}
