#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/FitSpline.h"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/MacroHelpers.h"
#include "../../includes/util/json_struct_def.hh"
#include <sstream>

using namespace ROOT;
using namespace ROOT::Experimental;

template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

static const char* label[] = {
	"21", "22", "23", "24", "41", "42", "31"
};
constexpr double zT = 3355 - 440/2;

using namespace hist;

enum class dir { x, y, both };

struct DoFit {
	struct Info { 
		int i; 
		dir d; 
		double lo, hi; 
	};
	std::vector<Info> data;
	DoFit() = default;
	DoFit(std::initializer_list<Info> lst) {
		for(auto elem : lst) data.push_back(elem);
	}
};

void tpc_draw_track (
	std::string fileName = "", 
	DoFit do_fit = {},
	A3 binning_x  = {100, -30, 30},
	A3 binning_xd = {100, -10, 10},
	A3 binning_y  = {100, -30, 30},
	A3 binning_yd = {100, -10, 10},
	A2 sci21_cut = {-DBL_MAX, DBL_MAX},
	A2 sci22_cut = {-DBL_MAX, DBL_MAX},
	A2 sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no
) {
	ROOT::EnableImplicitMT();
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_params = f->Get < 
			std::remove_reference_t<decltype(*tpc_params)>
		> ("FRS_tpc_parameters");
		if(!tpc_params)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
	}
	
	constexpr int N = 4;
	constexpr bool _take[N][2] = {
		{1, 1}, // 21
		{0, 0}, // 22
		{1, 1}, // 23
		{0, 0}  // 24
	};
	
	constexpr double WIDTH = 70.0;
	const Arr2<double, N, 2> zDL = [tpc_params](){
		Arr2<double,N,2> z{};
		for(int i=0; i<N; ++i) {
			z[i][0] = tpc_params->at(i).z0 - WIDTH/4;
			z[i][1] = tpc_params->at(i).z0 + WIDTH/4;
		}
		return z;
	}();
	const std::array<double, N> zTPC = {
		tpc_params->at(0).z0,	
		tpc_params->at(1).z0,	
		tpc_params->at(2).z0,	
		tpc_params->at(3).z0
	};
	printf("TPC positions: \n");
	for(int i=0; i<N; ++i) printf("TPC%s: %.1f mm\n", label[i], zTPC[i]);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60); 
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60); 
	auto* h2_xy = new TH2P(Form("Y [mm]:X [mm]@target z=%.1f", zT), 100, -40, 40, 100, -40, 40);
	auto* h2_ab = new TH2P(Form("Y-angle [mrad]:X-angle [mrad]@target z=%.1f", zT), 100, -20, 20, 100, -20, 20);

	TH2P* h2_tpc_xy[N];
	TH2P* h2_tpc_xd[N][2];
	TH2P* h2_tpc_yd[N][4];
	for(int i=0; i<N; ++i) {
		h2_tpc_xy[i] = new TH2P(Form("((h2_%d))Y measurement [mm]:X measurement [mm]@TPC%s", i, label[i]), 
			binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2]); 
		for(int dl: {0,1}) 
			h2_tpc_xd[i][dl] = new TH2P(Form("TPC%s X%d - Ref X [mm]:Ref X [mm]@TPC%s Delay line %d", label[i], dl, label[i], dl), 
				binning_x[0], binning_x[1], binning_x[2], binning_xd[0], binning_xd[1], binning_xd[2]);
		for(int a: {0,1,2,3})
			h2_tpc_yd[i][a] = new TH2P(Form("TPC%s Y%d - Ref Y [mm]:Ref Y [mm]@TPC%s Anode %d", label[i], a, label[i], a), 
				binning_y[0], binning_y[1], binning_y[2], binning_yd[0], binning_yd[1], binning_yd[2]);
	}
	
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	 
	std::vector<double> x; x.reserve(4*2);
	std::vector<double> y; y.reserve(4*2);
	std::vector<double> z; x.reserve(4*2);

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

		x.clear(); y.clear(); z.clear();
		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];
			const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;

			for(int d=0; d<2; ++d) {
				const std::vector<Measurement>& hits = tpc_hits[d];
				if(hits.size() == 1 and !std::isnan(hits[0].x) and _take[i][d]
					and !std::isnan(hits[0].y[0]) and !std::isnan(hits[0].y[1])) {
					x.push_back( hits[0].X() );
					y.push_back( hits[0].Y() );
					z.push_back( zDL[i][d] );
				}
			}

			h2_tpc_xy[i]->Fill( tpc.X0(), tpc.Y0() );
		}

		if(z.size() < 3) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	

		h2_xy->Fill(fx[1]*zT + fx[0], fy[1]*zT + fy[0]);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
	
		/* Now, for each TPC just fill the differences of individual measurements
		 * relative to the 'referent track' */
		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];
			
			for(int d=0; d<2; ++d) {
				const std::vector<Measurement>& hits = tpc.hits[d];
				if(hits.size() != 1) continue;
				
				const Measurement& m = hits[0];
				double x_ref = fx[1] * zDL[i][d] + fx[0];
				double y_ref = fy[1] * zDL[i][d] + fy[0];
				
				h2_tpc_xd[i][d] -> Fill( x_ref, m.x - x_ref );
				for(int a: {0,1}) {
					h2_tpc_yd[i][2*d + a] -> Fill( y_ref, m.y[a] - y_ref );
				}
			}
		}
	}

	TCanvas* cTr = new TCanvas("Tracks", "Tracks", 2000, 1200);
	cTr->Divide(2,2);
	cTr->cd(1);
	h2_track_x->Draw("COLZ");

	double r = 0.8;
	TLine* line;
	for(int i=0; i<4; ++i) {
		line = vline(h2_track_x, zTPC[i], r);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	line = vline(h2_track_x, zT, r);
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");

	cTr->cd(3);
	h2_track_y->Draw("COLZ");
	for(int i=0; i<4; ++i) {
		line = vline(h2_track_y, zTPC[i], r);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	line = vline(h2_track_y, zT, r);
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");
	cTr->cd(2); h2_xy->Draw("COLZ");
	cTr->cd(4); h2_ab->Draw("COLZ");
	
	/* ========================= */
	[[ maybe_unused ]] auto eval_latex_label = [&_take](int i) -> std::string {
		std::stringstream ref_txt {};
		if(_take[i][0] or _take[i][1]) {
			ref_txt << Form("TPC%s", label[i]);
			if(_take[i][0] and _take[i][1]) 
				ref_txt << " (both)";
			else if(!_take[i][0])
				ref_txt << " (DL2)";
			else
				ref_txt << " (DL1)";
		}
		return ref_txt.str();
	};

	[[ maybe_unused ]] std::string eval_code = [&_take]() -> std::string {
		std::stringstream ss;
		ss << "TPC_";
		for(int i=0; i<N; ++i) {
			int v = 0;
			v = static_cast<int>(_take[i][0]) |
				(static_cast<int>(_take[i][1]) << 1);
			if(v > 0) ss << label[i] << "(" << std::bitset<2>(v) << ")"; 
		}
		return ss.str();
	}();

	TCanvas* cdiff[N];
	for(int i=0; i<N; ++i) {
		TCanvas* c = new TCanvas(Form("TPC%s", label[i]), Form("TPC%s", label[i]), 1800, 1300); 
		c->Divide(3,2);
		for(int d: {0,1}) {
			c->cd(3*d + 1); 
			h2_tpc_xd[i][d]->Draw("COLZ");
		}
#define GET_ANODE_PAD(a) \
		{ \
			int d = a/2; \
			int ai = a%2; \
			c->cd( d*3 + 2 + ai ); \
		} 
		for(int a: {0,1,2,3}) {
			GET_ANODE_PAD(a);
			h2_tpc_yd[i][a]->Draw("COLZ");
		}
		//c->cd(4);
		//PLatex( 0.08, 
		//	Form("TPC%s", label[i]),
		//	"Reference made by: ",
		//	eval_latex_label(0),
		//	eval_latex_label(1),
		//	eval_latex_label(2),
		//	eval_latex_label(3)
		//);
		cdiff[i] = c;
	}
	TCanvas* cs = new TCanvas("SCIe", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	/* Idea is to take central bin positions (xi),
	 * their TH1D* projection has a peak around value (yi) with counts (wi)
	 * and then do a weighted linear fit. */
	const std::vector<DoFit::Info>& info = do_fit.data;

	for(const auto [i,o,lo,hi] : info) {
		TPCParam& p = tpc_params->at(i);
		TCanvas* c = cdiff[i];

		WARN("\n" EMPH1(TPC%s) " " EBOLD(fitting %s parameters) "\n", (o==dir::x?"X":(o==dir::y?"Y":"both")), label[i]);
		if(o == dir::x || o == dir::both) {
			for(int d: {0,1}) {
				const double b0 = p.x_offset[d];
				const double a0 = p.x_factor[d];

				auto [rg, gerr, g] = FitSpline<1, fit_info::GAUSS_MAX> (
						*h2_tpc_xd[i][d], lo, hi, 40, 1.1
						);
				auto [l,k] = rg;
				WARN("DL%d found \'more optimized\' parameter\n", d);
				WARN("\rOffset/slope: %.9f, %.9f\n", l, k);
				printf("\rBefore: (%.6f , %.6f)\n",             b0,                 a0);
				printf("\rNew   : " EBOLD((%.6f , %.6f)) "\n",  b0/(k+1) - l/(k+1), a0/(k+1));

				p.x_offset[d] = b0/(k+1) - l/(k+1);
				p.x_factor[d] = a0/(k+1);

				c->cd(3*d + 1); 
				gerr->Draw("P SAME");
				g->Draw("L SAME");
			}
		}
		if(o == dir::y || o == dir::both) {
			for(int a: {0,1,2,3}) {
				const double b0 = p.y_offset[a];
				const double a0 = p.y_factor[a];

				auto [rg, gerr, g] = FitSpline<1, fit_info::GAUSS_MAX> (
						*h2_tpc_yd[i][a], lo, hi, 40, 1.1
						);
				auto [l,k] = rg;
				WARN("A %d found \'more optimized\' parameter\n", a);
				WARN("\rOffset/slope: %.5f, %.5f\n", l, k);
				printf("\rBefore: (%.6f , %.6f)\n",             b0,                 a0);
				printf("\rNew   : " EBOLD((%.6f , %.6f)) "\n",  b0/(k+1) - l/(k+1), a0/(k+1));

				p.y_offset[a] = b0/(k+1) - l/(k+1);
				p.y_factor[a] = a0/(k+1);

				GET_ANODE_PAD(a); 
				gerr->Draw("P SAME");
				g->Draw("L SAME");
			}
		}
	}

	std::cout << setprecision(10);
	for(int i=0; i<4; ++i) {
		printf(BOLD ">> TPC%s: <<\n" KNRM, label[i]);
		std::cout << "\"x_factor\": " << tpc_params->at(i).x_factor << ",\n";
		std::cout << "\"x_offset\": " << tpc_params->at(i).x_offset << ",\n";
		std::cout << "\"y_factor\": " << tpc_params->at(i).y_factor << ",\n";
		std::cout << "\"y_offset\": " << tpc_params->at(i).y_offset << ",\n";
		std::cout << " ====================================================== \n";
	}

#if 0
	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { eval_code, inf.stem().c_str() });
	}
#endif
}
