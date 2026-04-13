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

enum class DoDiff {no, yes};
enum class DoOrientation {no, yes};
enum class DoSave {no, yes};

struct DoDelta {
	struct No {};
	struct Yes { 
		int position = -1;
		char orientation;
		double xlo, xhi;

		Yes() = default;
		constexpr Yes(int v) : position(v), orientation('\0'), xlo(0), xhi(0) {}; 
		constexpr Yes(int v, const char* o, double lo, double hi) : 
			position(v), orientation(o[0]), xlo(lo), xhi(hi) {}; 
	};

	/* constexpr */ static inline No no {};
	/* constexpr */ static inline Yes yes{-1};

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
	std::array<double,2> foot_cut = {10,3000},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	DoDelta do_delta = DoDelta::no,
	DoDiff do_diff = DoDiff::no,
	DoOrientation do_orientation = DoOrientation::no,
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

	const bool wrong_way = ((foot_param->orientation == "-x") or (foot_param->orientation == "-y"));
	const double offset = foot_param->delta_a;

	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};
	std::array<double, 4> zFOOT = {
		box->GetFOOTZ(0),
		box->GetFOOTZ(2),
		box->GetFOOTZ(4),
		box->GetFOOTZ(6)
	};

	TH2P* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	TH2P* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	TH2P* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);

	/* 4 possible placements along x, and y. */
	TH2P *foot_x[4];
	TH2P *foot_y[4];
	for(int i=0; i<4; ++i) {

		if(do_diff == DoDiff::no) {
			foot_x[i] = new TH2P(Form("((h2_x))FOOT Position Uncalibrated [mm]:TPC projection X at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]);
			foot_y[i] = new TH2P(Form("((h2_y))FOOT Position Uncalibrated [mm]:TPC projection Y at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]); 
		} else {
			foot_x[i] = new TH2P(Form("((h2_x))FOOT Position - TPC proj X [mm]:TPC projection X at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]);
			foot_y[i] = new TH2P(Form("((h2_y))FOOT Position - TPC proj Y [mm]:TPC projection Y at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]); 
		}
	}
	auto* foot_e_vs_d = new TH2P(Form("Cluster E [ADC Units G.M]:Delta [-0.5,0.5]@FOOT%d", ifoot), 
		80, -0.5, 0.5, 100, foot_cut[0], foot_cut[1]);
	auto* foot_e_vs_x = new TH2P(Form("Cluster E [ADC Units G.M.]:FOOT x@FOOT%d", ifoot), 
		binning_x[0], binning_x[1], binning_x[2], 100, foot_cut[0], foot_cut[1]);
	auto* foot_pos = new TH1P(Form("((h1))FOOT measurement@FOOT%d", ifoot), ORGB{0xFF7C0A}, binning_x[0], binning_x[1], binning_x[2]);
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
	
	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> z; x.reserve(N);
	
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

		double x_extrapolated, y_extrapolated;

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

		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		for(const auto& hit : foot->fCl) {
			double e = hit.fCE;
			double d = hit.Delta();
			if(do_delta == DoDelta::yes)
				e /= foot_param->de.CorrectionFactor(d);

			if(!mnd::IsInside(e, foot_cut)) continue;

			double hit_position = (hit.fCX - 319.5) * 0.150; // readout index to mm scale
			if(do_orientation == DoOrientation::yes and wrong_way)
				hit_position = -hit_position;

			hit_position -= offset;

			foot_pos->Fill(hit_position);
			for(int i=0; i<4; ++i) {
				x_extrapolated = fx[0] + fx[1]*zFOOT[i]; 
				y_extrapolated = fy[0] + fy[1]*zFOOT[i]; 
				
				if(do_diff == DoDiff::yes) {
					foot_x[i]->Fill(x_extrapolated, hit_position - x_extrapolated);
					foot_y[i]->Fill(y_extrapolated, hit_position - y_extrapolated);
				} else {
					foot_x[i]->Fill(x_extrapolated, hit_position);
					foot_y[i]->Fill(y_extrapolated, hit_position);
				}
			}
			foot_e_vs_d->Fill(d, e);
			foot_e_vs_x->Fill(hit_position, e);
		}
	}

	TCanvas* cTr = new TCanvas("TPC-tracks", "TPC-tracks", 2000, 1200);
	cTr->Divide(2,2);
	cTr->cd(1); h2_track_x->Draw("COLZ");
	cTr->cd(3); h2_track_y->Draw("COLZ");
	cTr->cd(2); h2_ab->Draw("COLZ");
	cTr->cd(4);
	PLatex(0.08,
		Form("Tracks derived from TPC: %s", ( []()->std::string {
			std::string s;
			for(int i=0; i<N; ++i) 
				if(_take[i]) { s+=::tpc_moniker[i]; s+=", "; };
			return s.substr(0, s.size()-2);
		}().c_str()))
	);

	TCanvas* c = new TCanvas(
		Form("FOOTXY%s", (do_diff == DoDiff::yes  ? "_diff" : "")), 
		Form("FOOT%d Position", ifoot),2000,1000);
	c->Divide(4,2);
	for(int i=0; i<4; ++i) {
		c->cd(i+1);
		foot_x[i]->Draw("COLZ");
		c->cd(i+5);
		foot_y[i]->Draw("COLZ");
	}

	if(do_delta == DoDelta::yes and do_delta.as_yes()->position != -1 and do_diff == DoDiff::yes) {
		int pos = do_delta.as_yes()->position;
		char o = do_delta.as_yes()->orientation;
		if(pos > 3) ERROR("Position %d requested, >3. Not allowed\n", pos);
		if(o != 'x' and o != 'y') ERROR("Orientation \'%c\' requested, isn't x or y.\n", o); 
		
		TH2D* hist; 
		if(o == 'x') { hist = &foot_x[pos]->h; c->cd(pos+1); }
		if(o == 'y') { hist = &foot_y[pos]->h; c->cd(pos+5); } 
		double xlo = do_delta.as_yes()->xlo;
		double xhi = do_delta.as_yes()->xhi;
		auto [alph, gerr, g] = FitSplineAndGraph<1, fit_info::GAUSS_MAX>(hist, xlo, xhi); 
		g->Draw("L SAME");
		gerr->Draw("P SAME");
		WARN("DoDelta(..) result: (%.4f, %.4f) meaning:\n"
			"Offset: %.4f, slope: %.5f\n", alph[0], alph[1], -alph[0]/(1+alph[1]), 1.0/(1+alph[1]));
	}

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* efoot = new TCanvas("cf", "FOOTE", 1400, 800);
	efoot->Divide(2,2);
	efoot->cd(1);
	foot_e_vs_d->Draw("COLZ"); gPad->SetLogz();
	efoot->cd(2);
	PLatex(0.07,
		[](){ std::string s = "Extrapolation done from: ";
			for(int i=0;i<N;++i) if(_take[i]) s += Form("TPC%s, ", tpc_moniker[i]);
			s.erase(s.size() - 2); return s; }(),
		Form("FOOT%d comes from DE10: %d", ifoot, foot_param->de10_index_),
		Form("In setup file it is placed at z0 = %.1f (total: %.1f)\n", 
			box->det_pos[ifoot], box->GetFOOTZ(ifoot, foot_param)),
		({ std::stringstream ss{}; ss << "Possible placements: " << box->det_pos; ss.str(); }),
		Form("Preliminary orientation: \'%s\'", foot_param->orientation.c_str()),
		"... Check if this matches!"
	);
	efoot->cd(3);
	foot_e_vs_x->Draw("COLZ"); gPad->SetLogz();
	printf("Positions at FOOT taken from extrapolating TPC: ");
	for(int i=0; i<N; ++i) {
		if(_take[i]) printf("%s ", tpc_moniker[i]);
	}
	printf("track!\n\n");

	efoot->cd(4);
	foot_pos->Draw();

	if(do_save == DoSave::yes)
		save_all(canvas::Extension::png, { Form("FOOT%d", ifoot) });
}
