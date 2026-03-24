#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../../includes/util/PolyFitter.hxx"
#include "../../includes/util/Tracking.hxx"
#include "../../includes/util/PrettyHisto.hxx"

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

enum class DoDelta {no, yes};
enum class DoDiff {no, yes};
enum class DoOrientation {no, yes};

/* Here we place FOOT ifoot onto various positions z, oriented either x- or y- 
 * and see which one fits the picture as correlation. */
void foot_spread (
	std::string fileName = "",
	uint32_t ifoot = 0, 
	std::array<double,2> foot_cut = {10,3000},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	DoDelta do_delta = DoDelta::no,
	DoDiff do_diff = DoDiff::no,
	DoOrientation do_orientation = DoOrientation::no
) {
	if(ifoot > 7) throw std::runtime_error("Second argument `ifoot` can be only {0,..,7}.");
	
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

	bool wrong_way = ((foot_param->orientation == "-x") or (foot_param->orientation == "-y"));
	
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

	/* 4 possible placements along x, and y. */
	TH2P *foot_x[4];
	TH2P *foot_y[4];
	for(int i=0; i<4; ++i) {

		if(do_diff == DoDiff::no) {
			foot_x[i] = new TH2P(Form("((h2_x))FOOT Position Uncalibrated [mm]:TPC projection X at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				300,-30, 30, 200, -30, 30);
			foot_y[i] = new TH2P(Form("((h2_y))FOOT Position Uncalibrated [mm]:TPC projection Y at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				300,-30, 30, 200, -30, 30); 
		} else {
			foot_x[i] = new TH2P(Form("((h2_x))FOOT Position - TPC proj X [mm]:TPC projection X at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				300,-30, 30, 400, -10, 10);
			foot_y[i] = new TH2P(Form("((h2_y))FOOT Position - TPC proj Y [mm]:TPC projection Y at Z=%.1f,pos=%d@FOOT%d", 
				zFOOT[i], i, ifoot),
				300,-30, 30, 400, -10, 10); 
		}
	}
	auto* foot_e = new TH2P(Form("Cluster E [ADC Units]:Delta [-0.5,0.5]@FOOT%d", ifoot), 
		80, -0.5, 0.5, 100, foot_cut[0], foot_cut[1]);
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

		for(const auto& hit : foot->fCl) {
			double e = hit.fCE;
			double d = hit.Delta();
			if(do_delta == DoDelta::yes)
				e /= foot_param->de.CorrectionFactor(d);

			if(!mnd::IsInside(e, foot_cut)) continue;

			double hit_position = (hit.fCX - 319.5) * 0.150; // readout index to mm scale
			if(do_orientation == DoOrientation::yes and wrong_way)
				hit_position = -hit_position;

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
			foot_e->Fill(d,e);
		}
	}

	TCanvas* c = new TCanvas("c", Form("FOOT%d Position", ifoot),2000,1000);
	c->Divide(4,2);
	for(int i=0; i<4; ++i) {
		c->cd(i+1);
		foot_x[i]->Draw("COLZ");
		c->cd(i+5);
		foot_y[i]->Draw("COLZ");
	}
	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 1800, 800);

	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* efoot = new TCanvas("cf", "FOOTE", 1400, 800);
	efoot->Divide(1,2);
	efoot->cd(1);
	foot_e->Draw("COLZ"); gPad->SetLogz();
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
	printf("Positions at FOOT taken from extrapolating TPC: ");
	for(int i=0; i<N; ++i) {
		if(_take[i]) printf("%s ", tpc_moniker[i]);
	}
	printf("track!\n\n");
}
