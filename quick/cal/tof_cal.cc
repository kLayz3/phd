/* Calibrate the time-of-flight between: 
 * [0] : SCI21 -> SCI22 , to indicate "possible" velocity at S2 
 * [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity at S3. 
 *
 * No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run. 
 * All of these measurements must be done with 12C files. There are 4 of relevance. */

#include "util/CLI.h"

#include "util/GaussFitMax.hxx"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.hxx"
#include "util/PolyFitter.h"
#include "util/Geometry.h"
#include "util/FitDrawer.hxx"
#include "util/Tracking.h"
#include "common/MacroCommon.hxx"

#include "IonOptics.hxx"

#include "TApplication.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::geom;

struct Brho {
    static constexpr char SEP = ';';
    double s2_incoming() const noexcept {
        return data_[0];
    }
    double s2_outgoing() const noexcept {
        return data_[1];
    };
    std::array<double, 2> data_;
};

struct FileBrho {
    static constexpr char SEP = ':';
    std::string name;
    Brho brho;
};
std::istream& operator>>(std::istream& , Brho& );
std::ostream& operator<<(std::ostream& , const Brho& );
std::istream& operator>>(std::istream& , FileBrho& );
std::ostream& operator<<(std::ostream& , const FileBrho& );

int main(int argc, char* argv[]) {
    CLI::App app{"Calibrate the time-of-flight between:\n\
                  [0] : SCI21 -> SCI22 , to indicate \"possible\" velocity at S2\n\
                  [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity and A/Q at S3.\n\
                  No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run (due to thick target).\n\
                  All of these measurements must be done with 12C files. There are 3 or 4 files of relevance."};
    
    
	std::vector<FileBrho> f;
	std::vector<TPCRef> ref{};
	std::array<double,3> dt_cut22_31 = {1000, -100, 100};
	std::array<double,3> dt_cut21_22 = {1000, -100, 100};
    double sratio = GAUSS_FIT_SIDE_RATIO_DEFAULT;
    double niter = 2;
    unsigned short line_size = 4;
	auto save = canvas::Extension::nil;
    add_logged_option(app, "-f,--file", f, "Pass one or more file names and corresponding 2 brho's.")
		->delimiter(',')
        ->type_name("[NAME:BRHO1;BRHO2 , ...]");
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
	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
	add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

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
	TApplication rootApp("app", 0, 0);

    constexpr static u32 Q0 = 6;
    constexpr static u32 A0 = 12;

    const u32 SCI21_I = RNFRSCal::SCI21_I;
    const u32 SCI22_I = RNFRSCal::SCI22_I;
    const u32 SCI31_I = RNFRSCal::SCI31_I;
    const u32 SCI_S2_I[2] = { SCI21_I, SCI22_I };
    const size_t nfiles = f.size();
    
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
    std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
    {
		std::unique_ptr<TFile> fhandle = std::make_unique<TFile>(f.front().name.c_str(), "READ");
        get_obj(fhandle, tpc_params, "FRS_tpc_parameters");
        get_obj(fhandle, sci_params, "FRS_sci_parameters");
    }
    constexpr auto N_TPC = TPCParam::N_S2_TPC;
	const Arr2<double, N_TPC, 2> zDL = TFRSCalCont::z_s2_tpc_delay_lines(tpc_params); 

    struct HistCont {
        TH1P *h_22_31, *h_21_22, *h_theta;
        TH2P *h_track_x;
    };
    std::vector<HistCont> hist;

	std::vector<double> x0, y0; // fitting containers ToF [0]
	std::vector<double> x1, y1; // fitting containers ToF [1]
    u32 i_clb_pnt = 1;
	for(const auto& [fileName, brho] : f) {
        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
        const double beta_inc = phy::Beta(Q0, A0, brho.s2_incoming());
        const double beta_out = phy::Beta(Q0, A0, brho.s2_outgoing());

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
		    mnd::PrintProgress(bar, entryId, nentries, 500, mnd::dancer2, 0.25);

            const auto& sci21 = frs->sci[SCI21_I];
            const auto& sci22 = frs->sci[SCI22_I];
			const auto& sci31 = frs->sci[SCI31_I];

            /* All ToF measuring stations must be single hit. No exception. */
			if(sci22.hits.size() != 1 or sci31.hits.size() != 1 or sci21.hits.size() != 1) 
                continue;
		    
            /* We want to measure the S2 angle,.. important for 21-22 ToF. */
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
        }
        hist.push_back( {
            .h_22_31 = h1_dt22_31,
            .h_21_22 = h1_dt21_22,
            .h_theta = h1_s2_angle,
            .h_track_x = h2_track_x } );
        auto [res22_31, _ ] = GaussFitMax(*h1_dt22_31, sratio, niter);
        auto [res21_22, __] = GaussFitMax(*h1_dt21_22, sratio, niter);
        x0.push_back( res22_31[1] ); // gauss peak value
        y0.push_back( 1.0 / beta_out );
        
        /* For S21-S22 β use their average. Normally, average path should also depend on theta,
         * but since it's mostly straight, and correction is O(theta^2), ignore it. */
        x1.push_back( res21_22[1] );
        y1.push_back( 2.0 / (beta_inc + beta_out) );

        bar.mark_as_completed();
        ++i_clb_pnt;
    }
    assert(hist.size() == nfiles);

    std::vector<double> fit_result_22_23, fit_result_21_22;
    auto [gerr0, g0] = FitAndDraw(1, x0, y0, {}, fit_result_22_23);
    auto [gerr1, g1] = FitAndDraw(1, x1, y1, {}, fit_result_21_22);
    gerr0->SetMarkerColor(kBlue - 1); gerr0->SetTitle("ToF Sci22 -> Sci31");
    gerr0->GetXaxis()->SetTitle("Mean ToF [25 ns], with offset");
    gerr0->GetYaxis()->SetTitle("1.0 / #beta");
    gerr1->SetMarkerColor(kRed - 1);  gerr1->SetTitle("ToF Sci21 -> Sci22");
    gerr1->GetXaxis()->SetTitle("Mean ToF [25 ns], with offset");
    gerr1->GetYaxis()->SetTitle("1.0 / #beta");
    g0->SetLineColor(kBlue);
    g1->SetLineColor(kRed);

    WARN("If the formula is: " EMPH(s / (ΔT + Λ))
         " then: " EMPH1(s = a ; Λ = -b\n));
    WARN("ToF S22 - S31: " EBOLD(b = %.6f; a = %.6f\n), fit_result_22_23[0], fit_result_22_23[1]);
    WARN("ToF S21 - S22: " EBOLD(b = %.6f; a = %.6f\n), fit_result_21_22[0], fit_result_21_22[1]);
    TCanvas *c = new TCanvas("Fit", "Fit", 1400, 700);
    c->Divide(2,1);
    c->cd(1); gerr0->Draw("AP"); g0->Draw("L SAME"); gPad->SetGrid();
    c->cd(2); gerr1->Draw("AP"); g1->Draw("L SAME"); gPad->SetGrid();

    TCanvas *c_raw = new TCanvas("RawToF", "RawToF", 2050, 1400);
    c_raw->Divide(4, nfiles);
    for(size_t i=0; i<nfiles; ++i) {
        auto [h_22_31, h_21_22, h_theta, h2_track_x] = hist[i];
        c_raw->cd(4*i + 1); h_22_31->DrawAndFit(sratio, kGreen, line_size, niter);
        c_raw->cd(4*i + 2); h_21_22->DrawAndFit(sratio, kGreen, line_size, niter);
        c_raw->cd(4*i + 3); h_theta->Draw();
        c_raw->cd(4*i + 4); h2_track_x->Draw("COLZ"); gPad->SetLogz();

        double r = 0.72;
        TLine* line;
        for(int sci: {0,1}) {
            line = hist::vline(h2_track_x, sci_params->at(SCI_S2_I[sci]).z0, r);
            line = hist::vline(h2_track_x, sci_params->at(SCI_S2_I[sci]).z0, r);
            line->SetLineColor(kRed);
            line->SetLineStyle(2);
            line->SetLineWidth(4);
            line->Draw("SAME");
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

std::istream& operator>>(std::istream& in, Brho& brho) {
    return ::mnd::template operator>> <Brho::SEP>(in, brho.data_);
}
std::ostream& operator<<(std::ostream& os, const Brho& brho) {
    return os << brho.data_;
}
std::istream& operator>>(std::istream& in, FileBrho& f) {
    if(!std::getline(in, f.name, FileBrho::SEP)) {
        WARN("Parsing Brho into string part (file name) failed.\n");
        return in;
    }

    return in >> f.brho;
}
std::ostream& operator<<(std::ostream& os, const FileBrho& f) {
    return os << f.name << ": " << f.brho; 
}
