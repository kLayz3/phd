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

enum class DoOffset {no, yes};

void foot_pair_angle (
	std::string fileName = "",
	uint32_t ipair = 0, 
	std::array<double,2> acceptance_x = {-30,30},
	std::array<double,2> acceptance_y = {-30,30},
	std::array<double,2> foot_x_cut = {3000,4000},
	std::array<double,2> foot_y_cut = {3000,4000},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	const size_t Npts = 1000,
	DoOffset do_offset = DoOffset::no,
	DoSave do_save = DoSave::no
) {
	if(ifoot > 3)
		throw std::invalid_argument("Second argument `ipair` can be only {0,1,2,3}.");
	
	using Measurement = RNTPCCal::Measurement;
	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	FOOTParam *foot_param_x, *foot_param_y; 
	FOOTBoxParam *box;
	{
		FOOTParam *p1, *p2;
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
		p1 = f->Get<FOOTParam>(Form("FOOT%d_setup", 2*ipair));
		p2 = f->Get<FOOTParam>(Form("FOOT%d_setup", 2*ipair + 1));
		if(!p1 || !p2)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
		box = f->Get<FOOTBoxParam>("FOOT0_box");
		if(!box)
			throw std::runtime_error(Form("FOOT box param is nullptr. Fix it (line: %d).", __LINE__));
		char o1 = p1->orientation[0];
		char o2 = p2->orientation[0];
		if(o1 == o2)
			ERROR("Orientations of subsequent entries have to be different\n");
		foot_param_x = (o1 == 'x') : p1 ? p2;
		foot_param_y = (o1 == 'y') : p1 ? p2;
	}

	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};
	
	const double z0 = box->GetFOOTZ(foot_param);
	const bool wrong_way = ((foot_param->orientation == "-x") or (foot_param->orientation == "-y"));
	const double offset = foot_param->delta_a;
	
	AngleFitResult::Direction df;
	if((foot_param->orientation == "-x") || (foot_param->orientation == "x")) 
		df = AngleFitResult::Direction::X;
	else 
		df = AngleFitResult::Direction::Y;

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
		if(sci21.hits.size() != 1) continue;
		if(sci22.hits.size() != 1) continue;
		if(sci31.hits.size() != 1) continue;
		
		h1_sci21->Fill(sci21.E);
		h1_sci22->Fill(sci22.E);
		h1_sci31->Fill(sci31.E);
		
		if(!mnd::IsInside(sci21.E, sci21_cut)) continue;
		if(!mnd::IsInside(sci22.E, sci22_cut)) continue;
		if(!mnd::IsInside(sci31.E, sci31_cut)) continue;

		h1_sci21_cut->Fill(sci21.E);
		h1_sci22_cut->Fill(sci22.E);
		h1_sci31_cut->Fill(sci31.E);

		double xRef_  = NAN;
		double yRef_  = NAN;
		double xFOOT_ = NAN;

		/* In an event, only a SINGLE valid FOOT value must be found. */
		bool is_foot_event_valid = false;
		for(const auto& hit : foot->fCl) {
			double e = hit.fCE;
			double d = hit.Delta();
			e /= foot_param->de.CorrectionFactor(d);
			
			if(!mnd::IsInside(e, foot_cut)) continue;
			
			double hit_position = (hit.fCX - 319.5) * 0.150; // readout index to mm scale
			if(wrong_way) 
				hit_position = -hit_position;

			hit_position -= offset;
			
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

		/* Extrapolated positions at the FOOT: */
		xRef_ = fx[1]*z0 + fx[0];
		yRef_ = fy[1]*z0 + fy[0];

		if(!mnd::IsInside(xRef_, acceptance_x)) continue; 
		if(!mnd::IsInside(yRef_, acceptance_y)) continue; 
		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_xy->Fill(xRef_, yRef_);
		h1_foot->Fill(xFOOT_);

		/* Add this point to the design vector(s) */
		xRef.push_back(xRef_);
		yRef.push_back(yRef_);
		xFOOT.push_back(xFOOT_);
		
		if(xFOOT.size() == Npts) {
			double angle;
			if(do_offset == DoOffset::yes) {
				auto r = FitAngleOffset(xRef, yRef, xFOOT);
				angle = r.t.Angle(df);
			} else {
				auto r = FitAngle(xRef, yRef, xFOOT);
				angle = r.Angle(df);
			}
			static int cnt = 0;
			double deg = angle * 180 / M_PI;
			printf("#%2d, phi = %.3f°\n", ++cnt, deg); 
			
			measurements.push_back(deg);
			xRef.clear(); yRef.clear(); xFOOT.clear();
		}
	}
	auto result = ::mean_stddev(measurements);
	std::cout << "=================\n"
		<< "Avg deg: " << result << "°\n";

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	
	TCanvas* cTr = new TCanvas("XY", "TPC-tracks", 2000, 1200);
	cTr->Divide(2,2);
	cTr->cd(1); h2_track_x->Draw("COLZ");
	cTr->cd(3); h2_track_y->Draw("COLZ");
	cTr->cd(2); h2_xy->Draw("COLZ");
	cTr->cd(4); h1_foot->Draw();

	TCanvas* cInfo = new TCanvas("info", "Info", 1600, 1000);
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

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("FOOT%d", ifoot) });
	}
}
