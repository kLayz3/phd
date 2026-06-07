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

constexpr const char* tpc_moniker[] = { "21", "22", "23", "24" };

constexpr int N = 4;
constexpr bool _take[N] = {
	1, // 21
	0, // 22
	1, // 23
	0  // 24
};

enum class DoOffset { no, yes};

void foot_angle (
	std::string fileName = "",
	uint32_t ifoot = 0, 
	std::array<double,3> binning_x = {200,-30,30},
	std::array<double,3> binning_y = {200,-30,30},
	std::array<double,2> foot_q_cut = {5.4, 6.6},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> angle_cut_x = {NAN, NAN},
	std::array<double,2> angle_cut_y = {NAN, NAN},
	const size_t Npts = 1000,
	DoOffset do_offset = DoOffset::no,
	DoSave do_save = DoSave::no
) {
	if(ifoot > 7)
		throw std::invalid_argument("Second argument `ifoot` can be only {0,..,7}.");
	
	using Measurement = RNTPCCal::Measurement;
	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	FOOTParam *foot_param; 
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
		foot_param = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!foot_param)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
		box = f->Get<FOOTBoxParam>("FOOT0_box");
		if(!box)
			throw std::runtime_error(Form("FOOT box param is nullptr. Fix it (line: %d).", __LINE__));
	}
	const Orientation o = foot_param->GetOrientation();
	if(o == Orientation::UNKNOWN)
		throw std::runtime_error("FOOT orientation not specified. I won't allow it.\n");

	const double z0 = box->GetFOOTZ(foot_param);
	const double R = foot_param->R();
	const char* ostr = (o == Orientation::X) ? "X" : "Y";
	const double offset = foot_param->delta_p;

	AngleFitResult::Direction df;
	if(o == Orientation::X)
		df = AngleFitResult::Direction::X;
	else 
		df = AngleFitResult::Direction::Y;
	
	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};

	WARN("Heuristically identified:\n" BOLD
		">> FOOT%d\n"
		">> DE10: %d\n"
		">> Measuring: \'%s\'\n"
		">> Inputted offset along measurement axis: %.2f mm\n"
		">> Z0: %.2f\n" KNRM,
		foot_param->N,
		foot_param->de10_index_,
		foot_param->orientation.c_str(),
		foot_param->delta_p,
		z0);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	TH2P* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent y-position [mm]:Referent X-position [mm]@FOOT%d", ifoot), 100, -40, 40, 100, -40, 40);
	auto* h1_foot = new TH1P("((h1_foot)) FOOT measurement [mm]", ORGB{0xC500CB}, 300, -50, 50);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	/* Containers for TPC extrapolation. */
	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> z; x.reserve(N);

	/* Containers for linear design problem. */
	std::vector<double> xRef;
	std::vector<double> yRef;
	std::vector<double> xFOOT;
	std::vector<double> measurements;

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
		
		x.clear(); y.clear(); z.clear();
		for(int i = 0; i < N; ++i) {
			if(i >= N or !_take[i]) continue;
			const auto& tpc = frs->tpc[i];

			double x0 = tpc.X0();
			double y0 = tpc.Y0();
			
			if(!std::isnan(x0) and !std::isnan(y0)) {
				x.push_back(x0);
				y.push_back(y0);
				z.push_back(zTPC[i]);
			}
		}
		if(z.size() < 2) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);
		
		if(mnd::IsValid(angle_cut_x) and !mnd::IsInside(fx[1]*1000, angle_cut_x)) continue;
		if(mnd::IsValid(angle_cut_y) and !mnd::IsInside(fy[1]*1000, angle_cut_y)) continue;

		/* In an event, only a single valid FOOT cluster must be found. */
		bool is_foot_event_valid = false;
		for(const auto& hit : foot->fCl) {
			double q = foot_param->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			double hit_position = foot_param->X0(hit);  // Takes care of offset as well. 
			
			if( std::isfinite(xFOOT_) ) {
				/* Already found valid point in the event. */
				is_foot_event_valid = false; break;
			} else {
				/* Export it outside. */
				xFOOT_ = hit_position;
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
			double angle;
			AngleFitResult r;
			if(do_offset == DoOffset::yes) {	
				auto tmp = FitAngleOffset(xRef, yRef, xFOOT);
				r = tmp.t;
			} else {
				r = FitAngle(xRef, yRef, xFOOT);
			}
			angle = r.Angle(df);
			static int cnt = 0;
			double deg = angle * 180 / M_PI;
			printf("#%2d, phi = %.3f°\n", ++cnt, deg); 
			
			measurements.push_back(deg);
			xRef.clear(); yRef.clear(); xFOOT.clear();
		}
	}
	auto result = mnd::mean_var(measurements);
	std::cout << "=================\n"
		<< "Avg deg: " << result << "°\n";

	TCanvas* cs = new TCanvas("SCIs&TPCs", "SCI21,22,31; TPC Ref", 1800, 800);
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
		Form("Referent track derived from TPC: %s", ( []()->std::string {
			std::string s;
			for(int i=0; i<N; ++i) 
				if(_take[i]) { s+=::tpc_moniker[i]; s+=", "; };
			return s.substr(0, s.size()-2);
		}().c_str())),
		Form("Number of measurements: %d", (int)measurements.size()),
		Form("Points per measurement: %zu", Npts),
		Form("Result: (%s)#circ", result.lstring().c_str())
	);
	cInfo->cd(3); h2_xy->Draw("COLZ");
	cInfo->cd(4); h2_ab->Draw("COLZ");

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("FOOT%d", ifoot) });
	}
}
