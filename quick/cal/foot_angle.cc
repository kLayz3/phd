#include "monad/accessories/formatting.hxx"
#include "util/CLI.h"
#include "util/PrettyHisto.h"
#include "util/MacroHelpers.h"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"
#include <numeric>

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

struct foot_enc {
	std::shared_ptr<RNFOOTCal> cont;
	FOOTParam* p;
	double z;
	double theta;
	double dxy;
	Orientation o;
};
std::ostream& operator<<(std::ostream& , const foot_enc& );

int main(int argc, char* argv[]) {
	CLI::App app{"\
Calibrate non-referent FOOT detectors' angle+offsets, based on the 9C beam test file run.\n\
The first two pairs' fully calibrated measurement serve as a referent measurement to draw\n\
particle tracks to match with other detectors in the sequence. "};
	
	constexpr u32 N_REF_PAIRS = 2;
    constexpr u32 N_REF_DET = 2*N_REF_PAIRS;

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
	add_logged_option(app, "-x,--bins-x", binning_x, "Binning X")
		->delimiter(',');
	add_logged_option(app, "-y,--bins-y", binning_y, "Binning Y")
		->delimiter(',');
	add_logged_option(app, "-q,--foot-cut",foot_q_cut, "FOOT Q cut (charge)")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21", sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22", sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31", sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
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
	if(Npts < 10)
		ERROR("Supply at least 10 points (-N,--npts). Is currently %zu.", Npts);
	
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
			rfoot[i].theta = rfoot[i].p->delta_a * M_PI / 180.0;
			rfoot[i].dxy = rfoot[i].p->delta_p;
			
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

	/* Verify that the referent detectors are indeed XY and properly set-up.
	 * There *must* be two X- and two Y-measuring detectors. */
	if(std::accumulate(std::begin(rfoot), std::end(rfoot), 0, [](int s, const auto& ref) {
		return s + (int)(ref.o == Orientation::X); }) != 2)
		ERROR("Expected %u X-measuring FOOTs in the %u referent pairs?", N_REF_PAIRS, N_REF_PAIRS);
	if(std::accumulate(std::begin(rfoot), std::end(rfoot), 0, [](int s, const auto& ref) {
		return s + (int)(ref.o == Orientation::Y); }) != 2)
		ERROR("Expected %u Y-measuring FOOTs in the %u referent pairs?", N_REF_PAIRS, N_REF_PAIRS);

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
	); {
		nlohmann::json j;
		j["delta_p"] = foot.p->delta_p;
		j["delta_a"] = foot.p->delta_a;
		std::cerr << "Before:\n" << MND_RGB_COL(201,129,255)
			<< Form("\"FOOT%u\": ", foot.p->de10_index_)
			<< j.dump(4) << KNRM << std::endl;
	}

	/* The normal-line: (x0,ax,y0,ay) converts to the measurement u|v as:
	 * u = ( cosθx_0     ,  zu*cosθx_0     , sinθx_0     , zu*sinθx_0) × (x0,ax,y0,ay) - offset
	 * v = (-sinθy_0     , -zv*sinθy_0     , cosθy_0     , zv*cosθy_0) × (x0,ax,y0,ay) - offset
	 */
	Eigen::Matrix<double, N_REF_DET, 4> A;
	Eigen::Vector<double, N_REF_DET> c;
	for(u32 i=0; i<N_REF_DET; ++i) {
		double t = rfoot[i].theta;
		double z = rfoot[i].z;
		if(rfoot[i].o == Orientation::X)
			A.row(i) <<  cos(t),  z*cos(t), sin(t), z*sin(t);
		else
			A.row(i) << -sin(t), -z*sin(t), cos(t), z*cos(t);

		c(i) = rfoot[i].dxy;
	}
	const auto qr = A.colPivHouseholderQr();
	if(qr.rank() < 4)
		ERROR("QR decomposition did not produce rank() >= 4. The referent FOOT part was:\n%s",
			mnd::streamable(rfoot).c_str());

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", 0x009B2F_c, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent Y-position [mm]:Referent X-position [mm]@FOOT%d", ifoot),
		binning_x[0],binning_x[1],binning_x[2],binning_y[0],binning_y[1],binning_y[2]);
	auto* h1_foot = new TH1P(Form("((h1_foot)) FOOT%u measurement [mm]", ifoot), 0xC500CB_c,
		binning[0],binning[1],binning[2]);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS");
	for(u32 i=0; i<N_REF_DET; ++i)
		rfoot[i].cont = model->MakeField<RNFOOTCal>(Form("FOOT%u", i));
    foot.cont = model->MakeField<RNFOOTCal>(Form("FOOT%u", ifoot));

	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	const size_t nentries = ntuple->GetNEntries();
	
	if(nentries < 1000) ERROR("ROOT file with less than 1k entries?");

	/* Optimally Npts has to be at least 3 times less than nentries */
	if(4*Npts > nentries) {
		WARN("4x Npts > nentries ; (4x %zu > %zu). Will trim -N,--npts to %zu",
			Npts, nentries, nentries/4);
		Npts = nentries/4;
	}

	/* Containers for referent detector measurements. */
	Eigen::Vector<double, N_REF_DET> ref_m;

	/* Containers for the final linear design problem. */
	std::vector<double> referent_xs;
	std::vector<double> referent_ys;
	std::vector<double> actual_measurements;

	/* Sequence of solved design problems. Stash their result here. */
	std::vector<double> result_a;
	std::vector<double> result_p;
    std::vector<double> result_s; // scale factor

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

		/* In an event, only a single valid FOOT cluster in the
         * to-be-calibrated detector must be found. */
		double m = NAN;
		u32 foot_n_counted = 0;
		for(const auto& hit : foot.cont->fCl) {
			const double q = foot.p->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			m = foot.p->BarePosition(hit);
			++foot_n_counted;
		}
		if(foot_n_counted != 1) continue;

		/* Find the reference containers. */
        bool ref_is_good = true;
		for(u32 i=0; i < N_REF_DET; ++i) {
			const auto& ft = rfoot[i];
		    u32 ncounted = 0;

            for(const auto& hit : ft.cont->fCl) {
				const double q = ft.p->Q( hit );
				if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;

				const double u_or_v = ft.p->BarePosition(hit);
				ref_m(i) = u_or_v;
				++ncounted;
			}
            /* In this case, just flag the entire event as unreliable. */
            if(ncounted != 1) {
                ref_is_good = false;
                break;
            }
		}
		if(!ref_is_good) continue;

		/* Solve the least squares problem:
		 * Finding the least squares solution of Ax = b is equivalent to solving the normal equation (AᵀA)x = Aᵀb.
		 * where b is ref_m + c. */
		
		const Eigen::Vector<double, N_REF_DET> b = ref_m + c;
		const Eigen::Vector4d refline = qr.solve(b);
		const auto [x0,ax,y0,ay] = refline;

		/* Extrapolated TRUE positions at the selected FOOT: */
		const double xRef = x0 + ax * z0;
		const double yRef = y0 + ay * z0;

		if(!mnd::IsInside(xRef, binning_x)) continue;
		if(!mnd::IsInside(yRef, binning_y)) continue;
		mnd::FillTrack(*h2_track_x, std::array{x0, ax});
		mnd::FillTrack(*h2_track_y, std::array{y0, ay});

		h2_ab->Fill(ax*1000.0, ay*1000); // [mrad]
		h2_xy->Fill(xRef, yRef); // is already in [mm].
		h1_foot->Fill(m);

		/* Add this point to the design problem vector(s) */
		referent_xs.push_back(xRef);
		referent_ys.push_back(yRef);
		actual_measurements.push_back(m);
		
		/* Once required statistics is reached, solve the design problem */
		if(actual_measurements.size() == Npts) {
			double angle, offset, scale_factor;

			const AngleOffsetFitResult r = FitAngleOffset(
				referent_xs,
				referent_ys,
				actual_measurements
			);
			
			angle = r.t.Angle(df);
			offset = r.c;
            scale_factor = r.t.a*r.t.a + r.t.b*r.t.b;

			static int count = 0;
			double deg = angle * 180 / M_PI;
			WARN("#%2d, phi = %.3f° ; gamma = %.3f;"
				"   a=%.4f, b=%.4f, a^2+b^2 = %.4f;\n",
				++count, deg, offset, r.t.a, r.t.b, scale_factor);
			
			result_a.push_back(deg);
			result_p.push_back(offset);
            result_s.push_back(scale_factor);

			referent_xs.clear();
			referent_ys.clear();
			actual_measurements.clear();
		}
	} // for(auto entryId : *ntuple)

	auto final_a = mnd::mean_var(result_a);
	auto final_p = mnd::mean_var(result_p);
    auto final_s = mnd::mean_var(result_s);

	/* p is what the fit gives, but this is the offset, so detector is placed at -mean */
	final_p.mean *= -1;
	std::cout << "\n========== TO BE WRITTEN INTO PARAM FILE =======" BOLD "\n"
		<< "Avg deg: " << final_a << "°\n"
		<< "Avg off: " << final_p << " mm\n"
        << Form("Scale factor (%u): ", ifoot) << final_s << KNRM "\n";
    nlohmann::json j;
    j["delta_p"] = final_p.mean;
    j["delta_a"] = final_a.mean;
    std::cerr << MND_RGB_BOLD_COL(207,255,182)
		<< Form("\"FOOT%u\": ", foot.p->de10_index_) << j.dump(4) 
		<< KNRM << std::endl;

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
            Form("Number of measurements: %d", (int)result_a.size()),
            Form("Points per measurement: %zu", Npts),
            Form("final angle: (%s)#circ", final_a.lstring().c_str()),
            Form("final offset: (%s) mm", final_p.lstring().c_str())
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
    return os 
		<< "z: " << foot.z
		<< "p: " << (foot.p ? mnd::ftos(foot.p->delta_p) + "mm" : "NULL")
		<< "a: " << (foot.p ? mnd::ftos(foot.p->delta_a) + "º" : "NULL")
		<< ", (" << magic_enum::enum_name(foot.o) << ')'
		<< " ";
}
