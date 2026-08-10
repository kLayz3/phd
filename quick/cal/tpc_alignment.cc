#include "util/CLI.h"

#include "util/FitSpline.h"
#include "util/PolyFitter.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.hxx"
#include "util/Tracking.h"
#include "common/MacroCommon.hxx"

#include "TApplication.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

struct HistDrawer {
    TH2P *hist = nullptr;
    TGraphErrors *gerr = nullptr;
    TGraph *g = nullptr;
    void Draw() {
        gPad->SetLogz();
        if(hist) hist->Draw("COLZ");
        if(g) g->Draw("L SAME");
        if(gerr) gerr->Draw("P SAME");
    }
};
int main(int argc, char* argv[]) {
	CLI::App app{"Calibrate the TPC detectors, by:\n\
                  => Defining referent (--ref) \"drift velocities\" for specific delay-lines (x-direction).\n\
                  => Defining referent drift velocities for specific anodes (y-direction).\n\
                  => Constraining the exact distances of each of the TPC's.\n\
                  \n\
                  All of these measurements must be done with 12C files. There are 3 files of relevance.\n\
                  To calibrate (the offsets) of the referent TPC delay lines/anodes, check the `tpc_alignment` in scripts."};

    constexpr auto N_TPC = TPCParam::N_S2_TPC;

    std::vector<std::string> fileNames{};
	std::vector<TPCRef> ref{}; 
	u32 i_tpc = 0;
	A3 binning_x = {100, -30, 30};
	A3 binning_y = {100, -30, 30};
	A3 binning_d = {50, -5, 5};
	A2 fit_range_x = {-5, 5};
	A2 fit_range_y = {-5, 5};
    u32 niter = 2;
    double sratio = 1.4;
	A2 sci21_cut = {-NAN, NAN};
	A2 sci22_cut = {-NAN, NAN};
	A2 sci31_cut = {-NAN, NAN};
    bool dont_fit = false;
	auto save = canvas::Extension::nil;

    add_logged_option(app, "-f,--file", fileNames, "Pass one or more file names from the ToF calibration runs.")
        ->delimiter(',')
        ->check(CLI::ReadPermissions);
    add_logged_option(app, "-i,--tpc", i_tpc, "Scintillator index; 0 => SCI21, 1 => SCI22.")
        ->check(CLI::Range(0, (int)N_TPC - 1));
    add_logged_option<DisplayDefault::No>(app, "-r, --ref", ref, 
		"Select which TPC's (either with index: 0,1,2, or with a label: 21,22,23) make the reference. \
		Select by '0/1' which delay lines get included into the measurement. ")
		->type_name("[INT|LABEL:BOOL,BOOL;...]")
		->delimiter(';');

	add_logged_option(app, "-x,--binning-x", binning_x, "Binning X.")
		->delimiter(',');
	add_logged_option(app, "-y,--binning-y", binning_y, "Binning Y.")
		->delimiter(',');
	add_logged_option(app, "-d,--binning-d", binning_d, "For the correlation plot, bining of the y-axis (TPC - ref), instead of the actual TPC measurement.")
        ->delimiter(',');
	add_logged_option(app, "--fit-range-x", fit_range_x, "Fit range [mm], along reference X-axis, given by binning-x option.")
		->delimiter(',');
	add_logged_option(app, "--fit-range-y", fit_range_y, "Fit range [mm], along reference Y-axis, given by binning-y option.")
		->delimiter(',');
    add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
        ->check(CLI::PositiveNumber);
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
        ->check(CLI::PositiveNumber);
    add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
    add_logged_flag(app, "--no-fit", dont_fit, "Do not fit the correlation plot.");
    add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.empty()) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

	if(ref.size() < 2) 
		ERROR("At least two valid referent TPC's must be given.\n");
    for(const auto& tpc : ref) {
		if(!tpc) { // operator bool() 
			std::cerr << tpc << std::endl; 
			ERROR("TPC invalid. Must be 0,1,2 and at least one dl flagged as valid."); 
		}
	}

    using Measurement = RNTPCCal::Measurement;
	TApplication rootApp("app", 0, 0);
	
    const auto& label = TFRSCalCont::tpc_label;
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
    {
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileNames.front().c_str(), "READ");
        get_obj(f, tpc_params, "FRS_tpc_parameters");
    }

	auto& tpc_param = tpc_params->at(i_tpc);
	const Arr2<double, N_TPC, 2> zDL = TFRSCalCont::z_s2_tpc_delay_lines(tpc_params); 
    const std::array<double, N_TPC> zTPC = TFRSCalCont::z_s2_tpc(tpc_params);
    const std::array<double, 2>& z0 = zDL[i_tpc];

    WARN("All TPC positions: \n");
	for(u32 i=0; i<N_TPC; ++i) printf("TPC%s: %.1f mm\n", label[i], zTPC[i]);
    
    WARN("Calibrating TPC%s: nominal positition: %.1f, delay lines nominally at: ",
         label[i_tpc], zTPC[i_tpc]);
    std::cerr << z0 << std::endl;
    WARN("Reference tracks are constructed by: " BOLD);
    for(const auto& r : ref) {
        fprintf(stderr, "TPC%s: %d|%d -- ", label[r.n], (int)r.use[0], (int)r.use[1]);
    } fprintf(stderr, KNRM "\n");

    TH2P *h2_x[2], *h2_y[4];
	for(int d=0; d<2; ++d)
		h2_x[d] = new TH2P(Form("TPC%s-X%d - Ref X [mm]:Ref X [mm]@At TPC%s", label[i_tpc], d, label[i_tpc]),
			binning_x[0], binning_x[1], binning_x[2], binning_d[0], binning_d[1], binning_d[2]);

	for(int a=0; a<4; ++a)
		h2_y[a] = new TH2P(Form("TPC%s-Y%d - Ref Y [mm]:Ref Y [mm]@At TPC%s", label[i_tpc], a, label[i_tpc]),
			binning_y[0], binning_y[1], binning_y[2], binning_d[0], binning_d[1], binning_d[2]);
	auto* h2_anode_mask = new TH2P("Anode mask [0,1,2,3]:Delay Line Index [0,1]", 4,-0.5,3.5, 2,-0.5,1.5);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
    
    auto* h2_track_x = new TH2P("((h2_track_x))Track density (X) [mm]:Depth z [mm]@S2 area", 800, 0, RNFRSCal::S2_LENGTH, 800, -60, 60);
    auto* h2_track_y = new TH2P("((h2_track_x))Track density (Y) [mm]:Depth z [mm]@S2 area", 800, 0, RNFRSCal::S2_LENGTH, 800, -60, 60);
    auto* h2_ab = new TH2P(Form("((h2_ab))Y-angle [mrad]:X-angle [mrad]@At TPC%s, clb point", label[i_tpc]), 100, -20, 20, 100, -20, 20);
    auto* h2_xy = new TH2P(Form("((h2_xy))Referent y-position [mm]:Referent X-position [mm]@At TPC%s, clb point", label[i_tpc]),
        binning_x[0],binning_x[1],binning_x[2],binning_y[0],binning_y[1],binning_y[2]);

    u32 f_index = 1;
    for(const auto& fileName : fileNames) {
        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
        
        ProgressBar bar {
            option::BarWidth{55},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg("[%u] TPC%s Alignment (per event: %s)", f_index, label[i_tpc], fileName.c_str())},
                option::ForegroundColor{Color::cyan},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
        };
        const size_t nentries = ntuple->GetNEntries();

        /* Containers for ref TPC extrapolation. */
        std::vector<double> xe, ye, ze;
        std::array<double, 2> xRef, yRef;

        for(auto entryId : *ntuple) {
            ntuple->LoadEntry(entryId);
            mnd::PrintProgress(bar, entryId, nentries, 500, mnd::dancer0, 0.20);
            
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

            /* Only look at hits when *both* delay lines have a good match. */
            const auto& tpc_to_calibrate = frs->tpc[i_tpc];
            const auto& hits_vec = tpc_to_calibrate.hits; // std::array< std::vector<..>, 2 >
            if(hits_vec[0].size() != 1 or hits_vec[1].size() != 1)
                continue;

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
                
            /* Extrapolated positions at the TPC, both delay lines: */
            for(u32 d: {0,1}) {
                xRef[d] = fx[1] * z0[d] + fx[0];
                yRef[d] = fy[1] * z0[d] + fy[0];
            }
            for(u32 d: {0,1}) {
                const Measurement hit = hits_vec[d].front();
                h2_x[d]->Fill (
                    xRef[d],
                    hit.X() - xRef[d] 
                );
                for(int a: {0,1})
                if( std::isfinite(hit.y[a]) ) {
                    h2_y[2*d + a]->Fill (
                        yRef[d],
                        hit.y[a] - yRef[d]
                    );
                }
                h2_anode_mask->Fill(d, hit.AnodeMask());
            }

            FillTrack(*h2_track_x, fx);
            FillTrack(*h2_track_y, fy);
            h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);
            h2_xy->Fill(mnd::mean(xRef), mnd::mean(yRef));
        }
        bar.mark_as_completed();
        ++f_index;
    }
    constexpr int N_PTS_FOR_GRAPH = 30;
    
    HistDrawer h_buff_dl[2], h_buff_a[4];
    if(!dont_fit) {
        TPCParam new_param = tpc_param;
        for(u32 d: {0,1}) {
            auto [rg, gerr, g] = FitSpline<1, fit_info::GAUSS_MAX> (
                *h2_x[d], fit_range_x[0], fit_range_x[1], N_PTS_FOR_GRAPH, sratio, niter /*, Verbosity::CHATTY */
            );
            auto [offset, slope]  = rg;
            new_param.x_offset[d] = (tpc_param.x_offset[d] - offset) / (slope + 1);
            new_param.x_factor[d] =       tpc_param.x_factor[d]      / (slope + 1);
            h_buff_dl[d] = { .hist = h2_x[d], .gerr = gerr, .g = g };
        }
        for(u32 a: {0,1,2,3}) {
            auto [rg, gerr, g] = FitSpline<1, fit_info::GAUSS_MAX> (
                *h2_y[a], fit_range_y[0], fit_range_y[1], N_PTS_FOR_GRAPH, sratio, niter /*, Verbosity::CHATTY */
            );
            auto [offset, slope]  = rg;
            new_param.y_offset[a] = (tpc_param.y_offset[a] - offset) / (slope + 1);
            new_param.y_factor[a] =       tpc_param.y_factor[a]      / (slope + 1);
            h_buff_a[a] = { .hist = h2_y[a], .gerr = gerr, .g = g };
        }

#define PRINT_TPC_PARAM(pname, pobj, delim) \
        fprintf(stderr, "\"%s\": ", #pname); std::cerr << (pobj).pname << #delim << std::endl;

        WARN("For TPC%s \n", label[i_tpc]);
        WARN("Old values: \n");
        PRINT_TPC_PARAM(x_factor, (tpc_param), (,))
        PRINT_TPC_PARAM(x_offset, (tpc_param), (,))
        PRINT_TPC_PARAM(y_factor, (tpc_param), (,))
        PRINT_TPC_PARAM(y_offset, (tpc_param), ())
        WARN("\nNew values: \n" BOLD)
        PRINT_TPC_PARAM(x_factor, (new_param), (,))
        PRINT_TPC_PARAM(x_offset, (new_param), (,))
        PRINT_TPC_PARAM(y_factor, (new_param), (,))
        PRINT_TPC_PARAM(y_offset, (new_param), ())

        fprintf(stderr, KNRM);
        for(int d: {0,1}) {
            const double slope0 = tpc_param.x_factor[d]; const double offset0 = tpc_param.x_offset[d];
            const double slope1 = new_param.x_factor[d]; const double offset1 = new_param.x_offset[d];
            WARN("~~~~ [DL%d] This is a relative change of: " BOLD "%.3f%% and %.3f%%\n" KNRM, d,
                100*std::abs((slope1-slope0)/slope0), 100*std::abs((offset1 - offset0)/offset0));
        }        
        for(int a: {0,1,2,3}) {
            const double slope0 = tpc_param.y_factor[a]; const double offset0 = tpc_param.y_offset[a];
            const double slope1 = new_param.y_factor[a]; const double offset1 = new_param.y_offset[a];
            WARN("~~~~ [AN%d] This is a relative change of: " BOLD "%.3f%% and %.3f%%\n" KNRM, a,
                100*std::abs((slope1-slope0)/slope0), 100*std::abs((offset1 - offset0)/offset0));
        }
    } else {
        for(u32 d: {0,1}) {
            h_buff_dl[d] = { .hist = h2_x[d], .gerr = nullptr, .g = nullptr };
        }
        for(u32 a: {0,1,2,3}) {
            h_buff_a[a] = { .hist = h2_y[a], .gerr = nullptr, .g = nullptr };
        }
    }
    
    using hist::vline;
	TCanvas* cTr = new TCanvas("Tracks", "Tracks", 2000, 1200);
	cTr->Divide(2,2);
	cTr->cd(1); gPad->SetLogz();
	h2_track_x->Draw("COLZ");

	double r = 0.8;
	TLine* line;
	for(int i=0; i<4; ++i) {
		line = vline(h2_track_x, zTPC[i], r);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}

	cTr->cd(3); gPad->SetLogz();
	h2_track_y->Draw("COLZ");
	for(int i=0; i<4; ++i) {
		line = vline(h2_track_y, zTPC[i], r);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	cTr->cd(2); gPad->SetLogz(); h2_xy->Draw("COLZ");
	cTr->cd(4); gPad->SetLogz(); h2_ab->Draw("COLZ");

    TCanvas* c = new TCanvas("Correlations", "Correlations", 2000, 1200);
    c->Divide(3,2);
    for(int d : {0,1}) {
        c->cd(3*d + 1); h_buff_dl[d].Draw();
        for(int a : {0,1}) {
            c->cd(3*d + 2 + a); h_buff_a[2*d + a].Draw();
        }
    }

    TCanvas* cs = new TCanvas("SCIe", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

    std::time_t now = std::time(nullptr);
    std::tm* tm = std::localtime(&now);
    char buffer[28];
    std::strftime(buffer, sizeof buffer,
        "%Y-%m-%d_%H-%M-%S", tm
    );
	canvas::save_all<canvas::Exe>(save, { Form("tpc%s", label[i_tpc]), buffer });

	WARN("End-of-main");
	rootApp.Run(); return 0;
}

