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

void verify_alignment (
	std::string fileName = "",
	uint32_t i = 0, 
	std::array<double,3> binning_x = {200,-10,10},
	std::array<double,3> binning_y = {200,-10,10},
	std::array<double,2> foot_q_cut = {5.4, 6.6},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> angle_cut_x = {NAN, NAN},
	std::array<double,2> angle_cut_y = {NAN, NAN},
	DoSave do_save = DoSave::no
) {
	using std::abs;

	TClass* cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!cl or !cl->GetDataMember("_x")) 
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	if(i >= N_PAIRS)
		ERROR("Second argument `i` can be only {0,..,%u}.", N_PAIRS);

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	std::array<FOOTParam,2>* foot_param; 
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param) ERROR("TPC param is nullptr.");
		foot_param = f->Get<std::array<FOOTParam,2>>(Form("FOOT_%d_setup", i));
		if(!foot_param) ERROR("FOOT_%d_setup is nullptr", i);
		
		box = f->Get<FOOTBoxParam>("FOOT_box");
		if(!box) ERROR("FOOT box param is nullptr.");
	}

	const double zfx_ = box->GetFOOTZ(&foot_param->at(0));
	const double zfy_ = box->GetFOOTZ(&foot_param->at(1));
	const double z0 = (zfx_ + zfy_) / 2;
	WARN("FOOTs heuristically identified as: X:%d, Y:%d\n"
		"Distances: X: %.1f | Y: %.1f\n"
		"All together: z0 = %.2f\n", 
		foot_param->at(0).N, foot_param->at(1).N,
		zfx_, zfy_, z0);

	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0,
		tpc_param->at(3).z0
	};

	auto* resx = new TH1P(Form("((h1_fitrx))DiffX FOOT-Ref[mm]@FOOT%d X", i),
			kMagenta-9, binning_x[0], binning_x[1], binning_x[2]);
	auto* resy = new TH1P(Form("((h1_fitry))DiffY FOOT-Ref [mm]@FOOT%d Y", i),
		kYellow-9, binning_y[0], binning_y[1], binning_y[2]);
	auto* resxy = new TH1P(Form("((h1_fitrxy))Total diff [mm]@FOOT%d-XY", i),
		kCyan-9, 
		(int)(binning_x[0] + binning_y[0])/2,
		0, 
		(abs(binning_x[1]) + abs(binning_y[1]) + abs(binning_x[2]) + abs(binning_y[2])) / 4
	);

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60);
	auto* h2_ab = new TH2P("Y-angle [mrad]:X-angle [mrad]", 100, -20, 20, 100, -20, 20);
	auto* h2_xy = new TH2P(Form("Referent y-position [mm]:Referent X-position [mm]@Pair%d", i), 
		320, -50, 50, 320, -50, 50);

	auto* h2_foot = new TH2P(Form("FOOT Y measurement [mm]:FOOT X measurement [mm]@Pair%d", i), 
		320, -50, 50, 320, -50, 50);

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
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);
		
		if(mnd::IsValid(angle_cut_x) and !mnd::IsInside(fx[1]*1000, angle_cut_x)) continue;
		if(mnd::IsValid(angle_cut_y) and !mnd::IsInside(fy[1]*1000, angle_cut_y)) continue;

		FillTrack(*h2_track_x, fx);
		FillTrack(*h2_track_y, fy);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		double x0 = fx[0] + fx[1]*z0;
		double y0 = fy[0] + fy[1]*z0; 
		h2_xy->Fill(x0, y0);
		
		for(const RNFOOTTrack& t : foot->track) {
			if(t.n != N_PAIRS) continue; // get only tracks going thru all 4 layers. 

			double  Q = t._q[i];
			double xf = t._x[i];
			double yf = t._y[i];
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(Q, foot_q_cut)) continue;
	
			h2_foot->Fill(xf, yf);
			resx->Fill(xf - x0);
			resy->Fill(yf - y0);
			
			double r = std::sqrt( (xf-x0)*(xf-x0) +  (yf-y0)*(yf-y0) );
			resxy->Fill(r);
		}
	}
	
	TCanvas* cs = new TCanvas("SCIs&TPCs", "SCI21,22,31; TPC Ref", 2150, 1400);
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

	TCanvas* cRes = new TCanvas("FitResidue", "Fit residue", 2150, 1400);
	cRes->Divide(3, 2);
	cRes->cd(1); resx->Draw();
	cRes->cd(2); resy->Draw();
	cRes->cd(3); resxy->Draw();
	cRes->cd(4); h2_xy->Draw("COLZ");
	cRes->cd(5); h2_foot->Draw("COLZ");
	cRes->cd(6);
	PLatex(0.08,
		"Here's some info..."
	);

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("pair%d",i) });
	}
}

