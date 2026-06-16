#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;

constexpr int N = 4;
constexpr bool _take[N] = {
	1, // 21
	0, // 22
	1, // 23
	0  // 24
};

enum class DoOffset { no, yes};

constexpr static A3 binning_da = {200, 0, 0.5};
void incoming_vs_outgoing (
	std::string fileName = "",
	A3 binning_x = {200,-10,10},
	A3 binning_y = {200,-10,10},
	A2 foot_q_cut = {5.4, 6.6},
	A2 sci21_cut = {-DBL_MAX, DBL_MAX},
	A2 sci22_cut = {-DBL_MAX, DBL_MAX},
	A2 sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no
) {
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param) ERROR("TPC param is nullptr.");
		
		box = f->Get<FOOTBoxParam>("FOOT_box");
		if(!box) ERROR("FOOT box param is nullptr.");
	}

	WARN(BOLD "Upstream reference constructed by: ");
	for(int i=0; i<N; ++i) if(_take[i]) fprintf(stderr, "TPC%s ", RNFRSCal::tpc_label[i]);
	fprintf(stderr, "\n");

	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};
	const double z0 = box->GetTargetZ();

	auto* diffx = new TH1P("Difference Outgoing - Incoming X [mm]",
		kMagenta-9, binning_x[0], binning_x[1], binning_x[2]);
	auto* diffy = new TH1P("Difference Outgoing - Incoming Y [mm]",
		kYellow-9, binning_y[0], binning_y[1], binning_y[2]);
	auto* diffxy = new TH1P("Total difference [mm]",
		kCyan-9, 
		(int)(binning_x[0] + binning_y[0])/2,
		0, 
		(abs(binning_x[1]) + abs(binning_y[1]) + abs(binning_x[2]) + abs(binning_y[2])) / 4
	);
	auto* diffa = new TH1P("Difference Angle [mrad]", kGreen-5, binning_da[0], binning_da[1], binning_da[2]);

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_ab_upst = new TH2P("TPC Y-angle [mrad]:TPC X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_ab_down = new TH2P("FOOT Y-angle [mrad]:FOOT X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy_upst = new TH2P("TPC y-position [mm]:TPC X-position [mm]@at target", 
		320, -50, 50, 320, -50, 50);
	auto* h2_xy_down = new TH2P("FOOT y-position [mm]:FOOT X-position [mm]@at target", 
		320, -50, 50, 320, -50, 50);
	auto* h2_diffa_vs_diffx = new TH2P("Angle diff [mrad]:Difference in X [mm]@FOOT - TPC at target",
		 binning_x[0], binning_x[1], binning_x[2], binning_da[0], binning_da[1], binning_da[2]);
	auto* h2_diffa_vs_diffy = new TH2P("Angle diff [mrad]:Difference in Y [mm]@FOOT - TPC at target",
		 binning_y[0], binning_y[1], binning_y[2], binning_da[0], binning_da[1], binning_da[2]);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSHit>("FRS"); 
	auto foot = model->MakeField<RNFOOTHit>("FOOT");
	auto ntuple = RNTupleReader::Open(std::move(model), "h104", fileName);

	/* Containers for TPC extrapolation. */
	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> z; x.reserve(N);

	ROOT::EnableImplicitMT();
	
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		const auto& sci21 = frs->cal.sci[0];
		const auto& sci22 = frs->cal.sci[1];
		const auto& sci31 = frs->cal.sci[2];
		
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
			if(i >= N or !_take[i]) continue;
			const auto& tpc = frs->cal.tpc[i];

			double x0 = tpc.X0();
			double y0 = tpc.Y0();
			
			if(!std::isnan(x0) and !std::isnan(y0)) {
				x.push_back(x0);
				y.push_back(y0);
				z.push_back(zTPC[i]);
			}
		}
		if(z.size() < 2) continue;

		if(foot->track.size() != 1) continue;
		const RNFOOTTrack& t = foot->track[0];
		if(t.n != N_PAIRS) continue; // get only tracks going thru all 4 layers.

		if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(t.Q, foot_q_cut)) continue;

		/* 3D line describing the track inside the box (downstream). */
		const mnd::geom::Line3D ft = RNTrackToLine3D(t);

		/* 2x 2D lines describing the track upstream. */
		const std::array<double,2> fx_tpc = PolyFit<1>(z, x);	
		const std::array<double,2> fy_tpc = PolyFit<1>(z, y);
		
		/* ======= TPC ======= */
		const double x0_tpc = fx_tpc[0] + fx_tpc[1]*z0;
		const double y0_tpc = fy_tpc[0] + fy_tpc[1]*z0; 
		h2_xy_upst->Fill(x0_tpc, y0_tpc);
		h2_ab_upst->Fill(fx_tpc[1]*1000, fy_tpc[1]*1000);
		
		FillTrack(*h2_track_x, fx_tpc, -HUGE_VAL, z0);
		FillTrack(*h2_track_y, fy_tpc, -HUGE_VAL, z0);
		
		/* ======= FOOT ======= */
		mnd::geom::Line2D fx_foot = ft.XLine();
		mnd::geom::Line2D fy_foot = ft.YLine();
		const double x0_foot = fx_foot.Eval( 0.0 );
		const double y0_foot = fy_foot.Eval( 0.0 );
		h2_xy_down->Fill(x0_foot, y0_foot);
		h2_ab_down->Fill(fx_foot[1]*1000, fy_foot[1]*1000);

		/* FOOT coordinate system is +z0 relative to the FRS one,
		 * represent it in the FRS coord. syst. which is `-z0` away from FOOTs' */
		fx_foot %= -z0;
		fy_foot %= -z0;
		FillTrack(*h2_track_x, fx_foot, z0, HUGE_VAL);
		FillTrack(*h2_track_y, fy_foot, z0, HUGE_VAL);

		double dx = x0_foot - x0_tpc; 
		double dy = y0_foot - y0_tpc; 
		double dr = std::sqrt( dx*dx + dy*dy );
		diffx ->Fill(dx);
		diffy ->Fill(dy);
		diffxy->Fill(dr);
		
		double da  = mnd::geom::Line3D(fx_foot, fy_foot)
			.AngleRelativeTo( mnd::geom::Line3D(fx_tpc,  fy_tpc) );

		diffa->Fill(da);
		h2_diffa_vs_diffx->Fill(dx, da);
		h2_diffa_vs_diffy->Fill(dy, da);
	}
	
	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2150, 1400);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* cTr = new TCanvas("Tracks", "Tracks", 2150, 1400);
	cTr->Divide(3,2);
	cTr->cd(1); h2_xy_upst->Draw("COLZ");
	cTr->cd(4); h2_xy_down->Draw("COLZ");
	cTr->cd(2); h2_ab_upst->Draw("COLZ");
	cTr->cd(5); h2_ab_down->Draw("COLZ");
	cTr->cd(3); h2_track_x->Draw("COLZ");
	cTr->cd(6); h2_track_y->Draw("COLZ");

	TCanvas* cDiff = new TCanvas("Diff", "FOOT-Reference differences", 2150, 1400);
	cDiff->Divide(3, 2);
	cDiff->cd(1); diffx->Draw();
	cDiff->cd(2); diffy->Draw();
	cDiff->cd(3); diffxy->Draw();
	cDiff->cd(4); diffa->Draw("COLZ");
	cDiff->cd(5); h2_diffa_vs_diffx->Draw("COLZ");
	cDiff->cd(6); h2_diffa_vs_diffy->Draw("COLZ");

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str() });
	}
}

