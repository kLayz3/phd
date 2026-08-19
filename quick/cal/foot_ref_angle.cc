#include "util/CLI.h"

#include "util/PrettyHisto.hxx"
#include "util/Tracking.h"
#include "util/MacroHelpers.h"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

#include "common/MacroCommon.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

struct foot_enc {
	std::shared_ptr<RNFOOTCal> cont;
	FOOTParam* p;
	double z;
	Orientation o;
};

int main(int argc, char* argv[]) {
	CLI::App app{"\
Calibrate referent FOOT detectors' (FOOT 0/1/2/3) angle, based on either the 12C or 9C beam calibration run.\n\
Due to the thick target, the position offsets shouldn't be touched directly, but rather make them 'agree'\n\
to referent 0 by manually fitting on the calibration run. This is done by a different program."};
	
	std::string fileName = "";
	u32 ifoot = 0;
	std::vector<TPCRef> ref{}; 
	std::array<double,3> binning_x = {200,-30,30};
	std::array<double,3> binning_y = {200,-30,30};
	std::array<double,2> foot_q_cut = {5.4, 6.6};
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX};
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX};
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX};
	size_t Npts = 1000;
	auto save = canvas::Extension::nil;

	add_logged_option(app, "-f,--file", fileName, "Pass a file name.")
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-i,--foot-id", ifoot, 
		"Select which FOOT detector.")
		->check(CLI::Range(0,3))
		->mandatory();
	add_logged_option<DisplayDefault::No>(app, "-r, --ref", ref, 
		"Select which TPC's (either with index: 0,1,2, or with a label: 21,22,23) make the reference. \
		Select by '0/1' which delay lines get included into the measurement. ")
		->type_name("[INT|LABEL:BOOL,BOOL;...]")
		->delimiter(';');
	add_logged_option(app, "-x,--bins-x",binning_x, "Binning X")
		->delimiter(',');
	add_logged_option(app, "-y,--bins-y",binning_y, "Binning Y")
		->delimiter(',');
	add_logged_option(app, "-q,--foot-cut",foot_q_cut, "FOOT Q cut (charge)")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option(app, "-N,--npts", Npts, "Amount of statistics required to issue a fit.");
	add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;
	if(ref.size() < 2) 
		ERROR("At least two valid referent TPC's must be given.\n");
	if(fileName.length() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}
	for(const auto& tpc : ref) {
		if(!tpc || !tpc.IsUpstream()) {
			ERROR("TPC%s (n=%u) flagged invalid. Must be upstream and at least one dl flagged as valid.",
                (tpc.n < RNFRSCal::N_VALID_TPC)? RNFRSCal::tpc_label[tpc.n]: "??", tpc.n); 
		}
	}

	TApplication rootApp("app", 0, 0);

	foot_enc foot; 
	FOOTBoxParam *box;
    std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");

        get_obj(f, box, "FOOT0_box");
		get_obj(f, foot.p,  Form("FOOT%u_setup", ifoot));
        get_obj(f, tpc_params, "FRS_tpc_parameters");
        
        // FOOT z assigned relative to FRS nominal Z, not the "box coordinates".
		foot.z = box->GetFOOTZ( foot.p );
		Orientation o = foot.p->GetOrientation();
		if(o == Orientation::UNKNOWN) ERROR("FOOT%d orientation not specified. I won't allow it.\n", ifoot);
		foot.o = o;
	}

    constexpr auto N_TPC = TPCParam::N_S2_TPC;
	const Arr2<double, N_TPC, 2> zDL = TFRSCalCont::z_s2_tpc_delay_lines(tpc_params); 
    const std::array<double, N_TPC> zTPC = TFRSCalCont::z_s2_tpc(tpc_params);

	const auto& binning = (foot.o == Orientation::X) ? binning_x : binning_y; 

	/* Reason these two enums are separate, 
	 * is so that the fitting library doesn't tag along any of the ROOT dependency. */ 
	auto df = (foot.o == Orientation::X)
		? AngleFitResult::Direction::X
		: AngleFitResult::Direction::Y;
		
	WARN("Heuristically identified:\n" BOLD
		">> FOOT%d\n" 
		">> DE10: %d\n"
		">> Measuring: \'%s\'\n"
		">> Inputted offset along measurement axis: %.2f mm\n"
		">> Zs: %.2f\n" KNRM,
		foot.p->N,
		foot.p->de10_index_,
		foot.p->orientation.c_str(),
		foot.p->delta_p,
		foot.z
	);

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
    auto* h2_track_x = new TH2P("((h2_track_x))Track density (X) [mm]:Depth z [mm]@S2 area", 800, 0, RNFRSCal::S2_LENGTH, 800, -60, 60);
    auto* h2_track_y = new TH2P("((h2_track_x))Track density (Y) [mm]:Depth z [mm]@S2 area", 800, 0, RNFRSCal::S2_LENGTH, 800, -60, 60);
	auto* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent y-position [mm]:Referent X-position [mm]@FOOT%d", ifoot), 
		binning_x[0],binning_x[1],binning_x[2],binning_y[0],binning_y[1],binning_y[2]);
	auto* h1_foot = new TH1P("((h1_foot)) FOOT measurement [mm]", ORGB{0xC500CB}, 
		binning[0],binning[1],binning[2]);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS");
	foot.cont = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	/* Containers for TPC extrapolation. */
	std::vector<double> xe, ye, ze;

	/* Containers for linear design problem. */
	std::vector<double> xRef;
	std::vector<double> yRef;
	std::vector<double> xFOOT;

	std::vector<double> measurements_a;
	std::vector<double> measurements_p;

	ROOT::EnableImplicitMT();
	
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

		double xRef_  = NAN;
		double yRef_  = NAN;
		double xFOOT_ = NAN;
		xe.clear(); ye.clear(); ze.clear();
		
		/* Find the reference containers. */
		for(const auto& id : ref) {
            u32 i = id.n;
			const auto& tpc = frs->tpc[i];
            for(u32 d : {0,1}) {
                if(!id.use[d] or tpc.hits[d].size() != 1)
                    continue;
                const RNTPCCal::Measurement& hits = tpc.hits[d].front();
                const double x = hits.X(); 
                const double y = hits.Y(); 
                if(!std::isfinite(x) or !std::isfinite(y)) 
                    continue;
                xe.push_back(x);
                ye.push_back(y);
                ze.push_back(zDL[i][d]);
            }
		}
		if(xe.size() < 3 or ye.size() < 3) continue;
		auto fx = PolyFit<1>(ze, xe);	
		auto fy = PolyFit<1>(ze, ye);

		/* In an event, only a single valid FOOT cluster must be found. */
		bool is_foot_event_valid = false;
		for(const auto& hit : foot.cont->fCl) {
			const double q = foot.p->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			const double hit_position = foot.p->BarePosition(hit); 
			
			if( std::isfinite(xFOOT_) ) {
				is_foot_event_valid = false; break; // Already found valid point in the event.
			} else {
				xFOOT_ = hit_position; // Export it outside.
			}
			is_foot_event_valid = true;
		}
		if(!is_foot_event_valid) continue;	

		/* Extrapolated positions at the FOOT: */
		xRef_ = fx[1] * foot.z + fx[0];
		yRef_ = fy[1] * foot.z + fy[0];

		if(!mnd::IsInside(xRef_, binning_x)) continue; 
		if(!mnd::IsInside(yRef_, binning_y)) continue; 
		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);
		h2_xy->Fill(xRef_, yRef_);
		h1_foot->Fill(xFOOT_);

		/* Add this point to the design problem vector(s) */
		xRef.push_back(xRef_);
		yRef.push_back(yRef_);
		xFOOT.push_back(xFOOT_);
		
		/* Once required statistics is reached, solve the design problem */
		if(xFOOT.size() == Npts) {
			double angle, offset;

			AngleOffsetFitResult r = FitAngleOffset(xRef, yRef, xFOOT);
			
			angle = r.t.Angle(df);
			offset = r.c;

			static int cnt = 0;
			double deg = angle * 180 / M_PI;
			WARN("#%2d, phi = %.3f° ; gamma = %.3f;"
				"   a=%.4f, b=%.4f, a^2+b^2 = %.4f;\n", 
				++cnt, deg, offset, r.t.a, r.t.b, r.t.a*r.t.a + r.t.b*r.t.b); 
			
			measurements_a.push_back(deg);
			measurements_p.push_back(offset);

			xRef.clear(); yRef.clear(); xFOOT.clear();
		}
	} // for(auto entryId : *ntuple)

	auto result_a = mnd::mean_var(measurements_a);
	auto result_p = mnd::mean_var(measurements_p);

	/* p is what the fit gives, but this is the offset, so detector is placed at -mean */
	result_p.mean *= -1;
	std::cout << "=================\n"
		<< BOLD "Avg deg: " << result_a << "°\n" << KNRM;

	TCanvas* cs = new TCanvas("SCIs&Refs", "SCI21,22,31; FOOT Ref", 2150, 1400);
	cs->Divide(3,3);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
    const double r = 0.8;
	TLine* line;
#define DRAW_LINES_TPC_AND_FOOT(hname) \
    { \
        for(int i=0; i < (int)TPCParam::N_S2_TPC; ++i) { \
            line = hist::vline(h2_track_x, zTPC[i], r); \
            line->SetLineColor(kRed); \
            line->SetLineStyle(2); \
            line->SetLineWidth(3); \
            line->Draw("SAME"); \
        } \
        line = hist::vline(h2_track_x, foot.z, r); \
        line->SetLineColor(kGreen + 1); \
        line->SetLineStyle(2); \
        line->SetLineWidth(4); \
        line->Draw("SAME"); \
    }
	cs->cd(7); h2_track_x->Draw("COLZ");
    DRAW_LINES_TPC_AND_FOOT(h2_track_x)
	cs->cd(8); h2_track_y->Draw("COLZ");
    DRAW_LINES_TPC_AND_FOOT(h2_track_y)
	
	TCanvas* cInfo = new TCanvas("Info", Form("Info-FOOT%d angle", ifoot), 2000, 1200);
	cInfo->Divide(2,2);
	cInfo->cd(1); h1_foot->Draw();
	cInfo->cd(2);
	new PLatex(0.08,
		"Referent track derived from FOOTs: 01&23",
		Form("Number of measurements: %d", (int)measurements_a.size()),
		Form("Points per measurement: %zu", Npts),
		Form("Result angle: (%s)#circ", result_a.lstring().c_str())
	);
	cInfo->cd(3); h2_xy->Draw("COLZ");
	cInfo->cd(4); h2_ab->Draw("COLZ");

	canvas::save_all<canvas::Exe>(save, { fileName, Form("FOOT%d", ifoot) });
	
	WARN("End-of-main");
	rootApp.Run(); return 0;
}
