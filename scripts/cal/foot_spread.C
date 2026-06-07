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

struct DoDelta {
	struct No {};
	struct Yes { 
		double xlo, xhi;
	};

	/* constexpr */ static inline No no {};
	/* constexpr */ static inline Yes yes {};

	DoDelta(No) : data_(No{}) {}
	DoDelta(Yes y) : data_(y) {}
	friend bool operator==(const DoDelta& lhs, const DoDelta& rhs) {
		return ( lhs.data_.index() == rhs.data_.index() &&
			lhs.data_.index() != std::variant_npos
		);
	}
	const Yes* as_yes() const { return std::get_if<Yes>(&data_); }

private:
	std::variant<No, Yes> data_;
};

/* Here we place FOOT ifoot onto various positions z, oriented either x- or y- 
 * and see which one fits the picture as correlation. */
void foot_spread (
	std::string fileName = "",
	uint32_t ifoot = 0, 
	std::array<double,3> binning_x = {200,-30,30},
	std::array<double,3> binning_y = {200,-30,30},
	std::array<double,2> foot_q_cut = {5.4, 6.6},
	std::array<double,2> sci21_cut = {NAN, NAN},
	std::array<double,2> sci22_cut = {NAN, NAN},
	std::array<double,2> sci31_cut = {NAN, NAN},
	std::array<double,2> angle_cut_x = {NAN, NAN},
	std::array<double,2> angle_cut_y = {NAN, NAN},
	DoDelta do_delta = DoDelta::no,
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
	Orientation o = foot_param->GetOrientation();
	if(o == Orientation::UNKNOWN)
		throw std::runtime_error("FOOT orientation not specified. I won't allow it.\n");

	const double z0 = box->GetFOOTZ(foot_param);
	const double R = foot_param->R();
	const char* ostr = (o == Orientation::X) ? "X" : "Y";
	
	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};

	TH2P* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	TH2P* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	TH2P* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);

	/* 4 possible placements along x, and y. */
	TH2P* foot_diff = new TH2P(Form("FOOT measurement - TPC ref [mm]:TPC projection %s at Z=%.1f,pos=%d@FOOT%d", 
		ostr, z0, ifoot, ifoot), binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]);
 
	auto* foot_q_vs_d = new TH2P(Form("Cluster Charge:Delta [-0.5, 0.5]@FOOT%d", ifoot), 
		80, -0.5, 0.5, 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_q_vs_x = new TH2P(Form("Cluster Charge:FOOT measurement [mm]@FOOT%d", ifoot), 
		binning_x[0], binning_x[1], binning_x[2], 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_pos = new TH1P(Form("((h1))FOOT measurement [mm]@FOOT%d", ifoot), ORGB{0xFF7C0A}, binning_x[0], binning_x[1], binning_x[2]);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
	/* Containers for TPC extrapolation. */
	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> z; x.reserve(N);
	
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

		double ref_extrapolated;

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
		double xFOOT_ = NAN;
		double qFOOT_ = NAN;
		double dFOOT_ = NAN;
		for(const auto& hit : foot->fCl) {
			
			double q = foot_param->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			double d = hit.Delta();
			double cx = hit.fCX;

			double hit_position = foot_param->BarePosition(hit);

			if( std::isfinite(xFOOT_) ) {
				/* Already found valid point in the event. */
				is_foot_event_valid = false; break;
			} else {
				/* Export it outside. */
				xFOOT_ = hit_position;
				qFOOT_ = q;
				dFOOT_ = d;
			}
			is_foot_event_valid = true;
		}
		if(!is_foot_event_valid) continue;

		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);
		foot_pos->Fill(xFOOT_);
		
		ref_extrapolated = (o == Orientation::X) ? (fx[0] + fx[1]*z0) : (fy[0] + fy[1]*z0); 
		
		foot_diff->Fill(ref_extrapolated, xFOOT_ - ref_extrapolated);
		foot_q_vs_d->Fill(dFOOT_, qFOOT_);
		foot_q_vs_x->Fill(xFOOT_, qFOOT_);
	}

	TCanvas* c = new TCanvas("FOOTDiff", Form("FOOT%d position difference plot", ifoot), 2000,1000);
	foot_diff->Draw("COLZ");

	double final_offset = 0.0;
	if(do_delta == DoDelta::yes) {
		TH2D* hist; 
		double xlo = do_delta.as_yes()->xlo;
		double xhi = do_delta.as_yes()->xhi;
		auto [alph, gerr, g] = FitSplineAndGraph<1, fit_info::GAUSS_MAX>(*foot_diff, xlo, xhi); 
		g->Draw("L SAME");
		gerr->Draw("P SAME");
		final_offset = -alph[0];
		WARN("DoDelta(..) result: (%.4f, %.4f) meaning:\n"
			"Offset to be put in JSON file: " EBOLD(%.4f) "\n", alph[0], alph[1], final_offset);
	}

	TCanvas* cs = new TCanvas("SCIs&TPCs", "SCI21,22,31 and TPC ref", 2200, 1200);
	cs->Divide(3,3);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	cs->cd(7); h2_track_x->Draw("COLZ");
	cs->cd(8); h2_track_y->Draw("COLZ");
	cs->cd(9); h2_ab->Draw("COLZ");

	TCanvas* efoot = new TCanvas("cf", "FOOTE", 1800, 800);
	efoot->Divide(2,2);
	efoot->cd(1);
	foot_q_vs_d->Draw("COLZ"); gPad->SetLogz();
	efoot->cd(2);
	PLatex(0.07,
		[](){ std::string s = "Extrapolation done from: ";
			for(int i=0;i<N;++i) if(_take[i]) s += Form("TPC%s, ", tpc_moniker[i]);
			s.erase(s.size() - 2); return s; }(),
		Form("FOOT%d comes from DE10: %d", ifoot, foot_param->de10_index_),
		Form("In setup file it is placed at z0 = %.1f (total: %.1f)\n", 
			box->det_pos[ifoot], z0),
		Form("Orientation: \'%s\', offset: %.4f mm", foot_param->orientation.c_str(),
			final_offset)
	);
	efoot->cd(3);
	foot_q_vs_x->Draw("COLZ"); gPad->SetLogz();
	printf("Positions at FOOT taken from extrapolating TPC: ");
	for(int i=0; i<N; ++i) {
		if(_take[i]) printf("%s ", tpc_moniker[i]);
	}
	printf("track!\n\n");

	efoot->cd(4);
	foot_pos->Draw();

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("FOOT%d", ifoot) });
	}
}
