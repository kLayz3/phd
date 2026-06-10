/* This script is meant to be run only after preliminary `foot_spread` has been
 * ran through. */

#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/FitSpline.hxx"
#include "../../includes/util/Tracking.hxx"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

struct foot_enc {
	std::shared_ptr<RNFOOTCal> cont;
	FOOTParam* p;
	double z;
	Orientation o;
};

void foot_angle (
	std::string fileName = "",
	uint32_t ifoot = 0, 
	std::array<double,3> binning_x = {200,-30,30},
	std::array<double,3> binning_y = {200,-30,30},
	std::array<double,2> foot_q_cut = {5.4, 6.6},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	const size_t Npts = 1000,
	DoSave do_save = DoSave::no
) {
	if(ifoot > 7 or ifoot < 4)
		throw std::invalid_argument("Second argument `ifoot` can be only {4,5,6,7}.");
	
	foot_enc foot[N_FOOT_DETECTORS]; 
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		box = f->Get<FOOTBoxParam>("FOOT0_box");
		if(!box) ERROR("FOOT box param is null");
		for(int i=0; i<N_FOOT_DETECTORS; ++i) {
			foot[i].p = f->Get<FOOTParam>(Form("FOOT%d_setup", i));
			if(!foot[i].p) ERROR("FOOT%d param is nullptr.", i);
		}
		double z_start = foot[0].z = box->GetFOOTZ( foot[0].p );
		for(int i=0; i<N_FOOT_DETECTORS; ++i) {
			foot[i].z = box->GetFOOTZ( foot[i].p ) - z_start;
			
			Orientation o = foot[i].p->GetOrientation();
			if(o == Orientation::UNKNOWN) ERROR("FOOT%d orientation not specified. I won't allow it.\n", i);
			foot[i].o = o;
		}
	}
	foot_enc& foot_to_calibrate = foot[ifoot];
	double z0 = foot_to_calibrate.z;
	const auto& binning = (foot_to_calibrate.o == Orientation::X) ? binning_x : binning_y; 

	/* Reason we keep it separate, is so that the fitting library doesn't tag along any of the ROOT dependency. */ 
	AngleFitResult::Direction df;
	if(foot_to_calibrate.o == Orientation::X)
		df = AngleFitResult::Direction::X;
	else 
		df = AngleFitResult::Direction::Y;
	
	WARN("Heuristically identified:\n" BOLD
		">> FOOT%d\n" 
		">> DE10: %d\n"
		">> Measuring: \'%s\'\n"
		">> Inputted offset along measurement axis: %.2f mm\n"
		">> Zs: ",
		foot_to_calibrate.p->N,
		foot_to_calibrate.p->de10_index_,
		foot_to_calibrate.p->orientation.c_str(),
		foot_to_calibrate.p->delta_p
	);
	for(int i=0; i<N_FOOT_DETECTORS; ++i) fprintf(stderr, "%.2f, ", foot[i].z);
	fprintf(stderr, KNRM "\n"); 

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, -20, box->width_outer+20, 600, -60, 60);
	auto* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent y-position [mm]:Referent X-position [mm]@FOOT%d", ifoot), 
		binning_x[0],binning_x[1],binning_x[2],binning_y[0],binning_y[1],binning_y[2]);
	auto* h1_foot = new TH1P("((h1_foot)) FOOT measurement [mm]", ORGB{0xC500CB}, 
		binning[0],binning[1],binning[2]);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS");
	for(int i=0; i<N_FOOT_DETECTORS; ++i)
		foot[i].cont = model->MakeField<RNFOOTCal>(Form("FOOT%d", i));
	
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	/* Containers for Pair0 & Pair1 extrapolation. */
	std::vector<double> xe, ye, zxe, zye;

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
		xe.clear(); ye.clear();	zxe.clear(); zye.clear();
		/* Find the reference containers. */
		for(int i=0; i<4; ++i) {
			const auto& ft = foot[i];
			for(const auto& hit : ft.cont->fCl) {
				double q = ft.p->Q( hit );
				if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;

				double hit_position = ft.p->X0(hit); // fully calibrated hit.
				if(ft.o == Orientation::X) {
					xe.push_back(hit_position);
					zxe.push_back( ft.z );
				} else {
					ye.push_back(hit_position);
					zye.push_back( ft.z );
				}
				break;
			}
		}
		if(xe.size() < 2 or ye.size() < 2) continue;
		auto fx = PolyFit<1>(zxe, xe);	
		auto fy = PolyFit<1>(zye, ye);
		//WARN("Fitting: "); std::cout << zxe << " : " << xe << " | " << zye << " : " << ye << std::endl;

		/* In an event, only a single valid FOOT cluster must be found. */
		bool is_foot_event_valid = false;
		for(const auto& hit : foot_to_calibrate.cont->fCl) {
			double q = foot_to_calibrate.p->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			double hit_position = foot_to_calibrate.p->BarePosition(hit); 
			
			if( std::isfinite(xFOOT_) ) {
				is_foot_event_valid = false; break; // Already found valid point in the event. 
			} else {
				xFOOT_ = hit_position; // Export it outside.
			}
			is_foot_event_valid = true;
		}
		if(!is_foot_event_valid) continue;	

		/* Extrapolated positions at the FOOT: */
		xRef_ = fx[1]*z0 + fx[0];
		yRef_ = fy[1]*z0 + fy[0];

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
	}
	auto result_a = mnd::mean_var(measurements_a);
	auto result_p = mnd::mean_var(measurements_p);

	/* p is what the fit gives, but this is the offset, so detector is placed at -mean */
	result_p.mean *= -1;
	std::cout << "=================\n"
		<< BOLD "Avg deg: " << result_a << "°\n"
		<< "Avg off: " << result_p << KNRM "\n";

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
	PLatex(0.08,
		"Referent track derived from FOOTs: 01&23",
		Form("Number of measurements: %d", (int)measurements_a.size()),
		Form("Points per measurement: %zu", Npts),
		Form("Result angle: (%s)#circ", result_a.lstring().c_str()),
		Form("Result offset: (%s) mm", result_p.lstring().c_str())
	);
	cInfo->cd(3); h2_xy->Draw("COLZ");
	cInfo->cd(4); h2_ab->Draw("COLZ");

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("FOOT%d", ifoot) });
	}
}
