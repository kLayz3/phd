#include "magic_enum/magic_enum.hpp"
#include "util/CLI.h"

#include "util/PrettyHisto.hxx"
#include "util/Tracking.h"
#include "util/MacroHelpers.h"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

struct foot_enc {
	std::shared_ptr<RNFOOTCal> cont;
	FOOTParam* p;
	double z;
	Orientation o;
};
std::ostream& operator<<(std::ostream& , const foot_enc& );

int main(int argc, char* argv[]) {
	CLI::App app{"\
Calibrate non-referent FOOT detectors' angle+offsets, based on the 9C beam test file run.\n\
The first two pairs' fully calibrated measurement serve as a referent measurement to draw\n\
particle tracks to match with other detectors in the sequence. "};
	
    constexpr static u32 N_REF_DET = 4;
	std::string fileName = "";
	u32 ifoot = 0;
	A3 binning_x = {200,-30,30};
	A3 binning_y = {200,-30,30};
	A2 foot_q_cut = {5.4, 6.6};
	A2 sci21_cut = {-DBL_MAX, DBL_MAX};
	A2 sci22_cut = {-DBL_MAX, DBL_MAX};
	A2 sci31_cut = {-DBL_MAX, DBL_MAX};
	size_t Npts = 1000;
    bool no_draw = false;
	auto save = canvas::Extension::nil;
    
	add_logged_option(app, "-f,--file", fileName, "Pass a file name.")
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-i,--foot-id", ifoot, 
		"Select which FOOT detector.")
		->check(CLI::Range((int)N_REF_DET, N_FOOT_DETECTORS-1))
		->mandatory();
	add_logged_option(app, "-x,--bins-x",binning_x, "Binning X")
		->delimiter(',');
	add_logged_option(app, "-y,--bins-y",binning_x, "Binning Y")
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
    add_logged_flag(app, "--nodraw", no_draw, "Do not draw the histograms.");
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileName.length() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}
	
	TApplication rootApp("app", 0, 0);

	foot_enc rfoot[N_REF_DET];
    foot_enc foot;
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
        get_obj(f, box, "FOOT0_box");
		box = f->Get<FOOTBoxParam>("FOOT0_box");
		if(!box) ERROR("FOOT box param is null");
		for(u32 i=0; i<N_REF_DET; ++i) {
            get_obj(f, rfoot[i].p, Form("FOOT%u_setup", i));
			rfoot[i].z = box->GetFOOTZRel( rfoot[i].p );
			Orientation o = rfoot[i].p->GetOrientation();
			if(o == Orientation::UNKNOWN) ERROR("Referent FOOT%u orientation not specified. I won't allow it.\n", i);
			rfoot[i].o = o;
		}
        get_obj(f, foot.p, Form("FOOT%u_setup", ifoot));
        foot.z = box->GetFOOTZRel( foot.p );
        Orientation o = foot.p->GetOrientation();
        if(o == Orientation::UNKNOWN) ERROR("To-calibrate FOOT%d orientation not specified. I won't allow it.\n", ifoot);
        foot.o = o;
	}
	double z0 = foot.z;
	const auto& binning = (foot.o == Orientation::X) ? binning_x : binning_y;

	/* Reason these two enums are separate,
	 * is so that the fitting library doesn't tag along any of the ROOT dependency. */
	auto df = (foot.o == Orientation::X)
		? AngleFitResult::Direction::X
		: AngleFitResult::Direction::Y;
	
    WARN("Referent FOOT detectors given as:\n" );
	for(u32 i=0; i<N_REF_DET; ++i) {
        std::cerr << rfoot[i] << "\n";
    }
	WARN("Wanted detector heuristically identified:\n" BOLD
		">> FOOT%d\n"
		">> DE10: %d\n"
		">> Measuring: \'%s\'\n"
		">> Inputted offset along measurement axis: %.2f mm\n"
		">> Z: %.4f\n" KNRM,
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
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent Y-position [mm]:Referent X-position [mm]@FOOT%d", ifoot),
		binning_x[0],binning_x[1],binning_x[2],binning_y[0],binning_y[1],binning_y[2]);
	auto* h1_foot = new TH1P(Form("((h1_foot)) FOOT%u measurement [mm]", ifoot), ORGB{0xC500CB},
		binning[0],binning[1],binning[2]);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS");
	for(u32 i=0; i<N_REF_DET; ++i)
		rfoot[i].cont = model->MakeField<RNFOOTCal>(Form("FOOT%u", i));
    foot.cont = model->MakeField<RNFOOTCal>(Form("FOOT%u", ifoot));

	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	/* Containers for referent detector, Pair0 & Pair1 extrapolation. */
	std::vector<double> xe, ye, zxe, zye;

	/* Containers for linear design problem. */
	std::vector<double> xRef;
	std::vector<double> yRef;
	std::vector<double> xFOOT;

	std::vector<double> measurements_a;
	std::vector<double> measurements_p;
    std::vector<double> measurements_s; // scale factor

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
		double mFOOT_ = NAN;
		xe.clear(); ye.clear();	zxe.clear(); zye.clear();
		
		/* Find the reference containers. */
        bool ref_is_good = true;
		for(u32 i=0; i < N_REF_DET; ++i) {
			const auto& ft = rfoot[i];
		    u32 ncounted = 0;
            
            for(const auto& hit : ft.cont->fCl) {
				const double q = ft.p->Q( hit );
				if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;

				double hit_position = ft.p->X0(hit); // fully calibrated hit.
				if(ft.o == Orientation::X) {
					xe.push_back(hit_position);
					zxe.push_back( ft.z );
				} else {
					ye.push_back(hit_position);
					zye.push_back( ft.z );
				}
				++ncounted;
			}
            /* In this case, just flag the entire event as unreliable. */
            if(ncounted != 1) {
                ref_is_good = false;
                break;
            }
		}
		if(!ref_is_good or xe.size() < 2 or ye.size() < 2) continue;

		/* In an event, only a single valid FOOT cluster in the
         * to-be-calibrated detector must be found. */
		u32 foot_n_counted = 0;
        mFOOT_ = NAN; // just in case
		for(const auto& hit : foot.cont->fCl) {
			const double q = foot.p->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			mFOOT_ = foot.p->BarePosition(hit);
			++foot_n_counted;
		}
		if(foot_n_counted != 1) continue;

		const auto fx = PolyFit<1>(zxe, xe);
		const auto fy = PolyFit<1>(zye, ye);
		/* Extrapolated positions at the selected FOOT: */
		xRef_ = fx[1]*z0 + fx[0];
		yRef_ = fy[1]*z0 + fy[0];

		if(!mnd::IsInside(xRef_, binning_x)) continue;
		if(!mnd::IsInside(yRef_, binning_y)) continue;
		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);
		h2_xy->Fill(xRef_, yRef_);
		h1_foot->Fill(mFOOT_);

		/* Add this point to the design problem vector(s) */
		xRef.push_back(xRef_);
		yRef.push_back(yRef_);
		xFOOT.push_back(mFOOT_);
		
		/* Once required statistics is reached, solve the design problem */
		if(xFOOT.size() == Npts) {
			double angle, offset, scale_factor;

			AngleOffsetFitResult r = FitAngleOffset(xRef, yRef, xFOOT);
			
			angle = r.t.Angle(df);
			offset = r.c;
            scale_factor = r.t.a*r.t.a + r.t.b*r.t.b;

			static int count = 0;
			double deg = angle * 180 / M_PI;
			WARN("#%2d, phi = %.3f° ; gamma = %.3f;"
				"   a=%.4f, b=%.4f, a^2+b^2 = %.4f;\n",
				++count, deg, offset, r.t.a, r.t.b, scale_factor);
			
			measurements_a.push_back(deg);
			measurements_p.push_back(offset);
            measurements_s.push_back(scale_factor);

			xRef.clear(); yRef.clear(); xFOOT.clear();
		}
	} // for(auto entryId : *ntuple)

	auto result_a = mnd::mean_var(measurements_a);
	auto result_p = mnd::mean_var(measurements_p);
    auto result_s = mnd::mean_var(measurements_s);

	/* p is what the fit gives, but this is the offset, so detector is placed at -mean */
	result_p.mean *= -1;
	std::cout << "========== TO BE WRITTEN INTO PARAM FILE =======" BOLD "\n"
		<< "Avg deg: " << result_a << "°\n"
		<< "Avg off: " << result_p << " mm\n"
        << Form("Scale factor (%u): ", ifoot) << result_s << KNRM "\n";
    nlohmann::json j;
    j["delta_p"] = result_p.mean;
    j["delta_a"] = result_a.mean;
    std::cerr << Form("\"FOOT%u\": ", foot.p->de10_index_) << j.dump(4) << std::endl;

    if(!no_draw) {
        TCanvas* cs = new TCanvas("SCIs&Refs", "SCI21,22,31; FOOT Ref", 2150, 1400);
        cs->Divide(3,3);
        cs->cd(1); h1_sci21->Draw();
        cs->cd(2); h1_sci22->Draw();
        cs->cd(3); h1_sci31->Draw();
        cs->cd(4); h1_sci21_cut->Draw();
        cs->cd(5); h1_sci22_cut->Draw();
        cs->cd(6); h1_sci31_cut->Draw();
        cs->cd(7); h2_track_x->Draw("COLZ");
        cs->cd(8); h2_track_y->Draw("COLZ");
        
        TCanvas* cInfo = new TCanvas("Info", Form("Info-FOOT%d angle", ifoot), 2000, 1200);
        cInfo->Divide(2,2);
        cInfo->cd(1); h1_foot->Draw();
        cInfo->cd(2);
        new PLatex(0.08,
            "Referent track derived from FOOTs: 01&23",
            Form("Number of measurements: %d", (int)measurements_a.size()),
            Form("Points per measurement: %zu", Npts),
            Form("Result angle: (%s)#circ", result_a.lstring().c_str()),
            Form("Result offset: (%s) mm", result_p.lstring().c_str())
        );
        cInfo->cd(3); h2_xy->Draw("COLZ");
        cInfo->cd(4); h2_ab->Draw("COLZ");

        canvas::save_all<canvas::Exe>(save, { Form("FOOT%d", ifoot) });
        
        WARN("End-of-main");
        rootApp.Run();
    }
    return 0;
}

std::ostream& operator<<(std::ostream& os, const foot_enc& foot) {
    return os << foot.z << '(' << magic_enum::enum_name(foot.o) << ')';
}
