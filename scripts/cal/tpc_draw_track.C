#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../../includes/util/PolyFitter.hxx"
#include "../../includes/util/Tracking.hxx"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/json_struct_def.hh"

using namespace ROOT;
using namespace ROOT::Experimental;

template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

static const char* label[] = {
	"21", "22", "23", "24", "41", "42", "31"
};
constexpr double zT = 3355 - 440/2;

double lo_y(TH2D*, double );
double hi_y(TH2D*, double );

std::array<double, 2> FitProfile(TH2D*, int);
std::array<double, 2> FitProfile(TH2D*, double);

struct OptionE {
	static constexpr int Nothing = 0x0000;
	static constexpr int FitDiff = 0x0001;
	int mask_;
	inline bool operator==(int rhs) {
		switch(rhs) {
			case 0:  return mask_ == rhs;
			default: return static_cast<bool>( mask_ & rhs );
		}
	}
	OptionE() = default;
	OptionE(int x) : mask_(x) {}
};

void tpc_draw_track(std::string fileName = "", OptionE e = OptionE::Nothing, int I = -1, uint64_t max_events = -1) {
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

	TH2P* h2_xy = new TH2P(Form("Y [mm]:X [mm]@target z=%.1f", zT), 400, -50, 50, 400, -50, 50);
	TH2P* h2_ab = new TH2P(Form("Y-angle [mrad]:X-angle [mrad]@target z=%.1f", zT), 400, -50, 50, 400, -50, 50);
	TH2P* h2_track_x = new TH2P("Track density (X) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60); 
	TH2P* h2_track_y = new TH2P("Track density (Y) [mm]:Depth z [mm]@S2 area", 600, 0, 4500, 500, -60, 60); 

	TH2P* h2_tpc_xy[N];
	TH2P* h2_tpc_xd[N][2];
	TH2P* h2_tpc_yd[N][4];
	for(int i=0; i<N; ++i) {
		h2_tpc_xy[i] = new TH2P(Form("TPC%s Y [mm]:TPC%s X [mm]", label[i], label[i]), 400, -20, 20, 400, -20, 20); 
		for(int dl: {0,1}) 
			h2_tpc_xd[i][dl] = new TH2P(Form("TPC%s X%d - Ref Y [mm]:Ref X [mm]@Delay line %d", label[i], dl, dl), 200,-50,50,200,-10,10);
		for(int a: {0,1,2,3})
			h2_tpc_yd[i][a] = new TH2P(Form("TPC%s Y%d - Ref Y [mm]:Ref Y [mm]@Anode %d", label[i], a, a), 200,-100,100,200,-10,10);
	}
	constexpr int N_BINS_X = 3000;
	constexpr int X_LO     = -30;
	constexpr int X_HI     =  30;
	constexpr int N_BINS_XD = 150;
	constexpr int XD_LO     = -10;
	constexpr int XD_HI     =  10;

	constexpr int N_BINS_Y = 300;
	constexpr int Y_LO     = -50;
	constexpr int Y_HI     =  50;
	constexpr int N_BINS_YD = 200;
	constexpr int YD_LO     = -10;
	constexpr int YD_HI     =  10;

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	 
	std::vector<double> x; x.reserve(4*2);
	std::vector<double> y; y.reserve(4*2);
	std::vector<double> z; x.reserve(4*2);
	
	uint64_t maxEntries = std::min( max_events, static_cast<uint64_t>(ntuple->GetNEntries()) );

	for(uint64_t entryId = 0; entryId < maxEntries; ++entryId) {
		ntuple->LoadEntry(entryId);
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
				double x_ref = fx[1]*zDL[i][d] + fx[0];
				double y_ref = fy[1]*zDL[i][d] + fy[0];
				
				h2_tpc_xd[i][d] -> Fill( x_ref, m.x - x_ref );
				for(int a: {0,1}) {
					h2_tpc_yd[i][2*d + a] -> Fill( y_ref, m.y[a] - y_ref );
				}
			}
		}
	}

	TCanvas* cT = new TCanvas("cT", "cT", 2000, 1200);
	cT->Divide(1,2);
	cT->cd(1);
	h2_xy->Draw("COLZ");

	cT->cd(2);
	h2_ab->Draw("COLZ");

	/* ========================= */
	/* Also draw some nominal indicators where each TPC sits. */

	TCanvas* cTr = new TCanvas("cTr", "cTr", 2000, 1200);
	cTr->Divide(1,2);
	cTr->cd(1);
	h2_track_x->Draw("COLZ");

	double r = 0.8;
	for(int i=0; i<4; ++i) {
		TLine* line = new TLine( zTPC[i], lo_y(*h2_track_x, r), zTPC[i], hi_y(*h2_track_x, r) );
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	TLine* line = new TLine( zT, lo_y(*h2_track_x, r), zT, hi_y(*h2_track_x, r) );
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");

	cTr->cd(2);
	h2_track_y->Draw("COLZ");
	for(int i=0; i<4; ++i) {
		TLine* line = new TLine( zTPC[i], lo_y(*h2_track_y, r), zTPC[i], hi_y(*h2_track_y, r) );
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	line = new TLine( zT, lo_y(*h2_track_y, r), zT, hi_y(*h2_track_y, r) );
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");

	/* ========================= */
	std::stringstream ref_txt{};
	auto eval_latex_label = [&_take](int i) -> std::string {
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

	for(int i=0; i<N; ++i) {
		if(I >= 0 and i != I) continue;

		TCanvas* c = new TCanvas(Form("TPC%s", label[i]), Form("TPC%s", label[i]), 1800, 1300); 
		c->Divide(4,2);
		for(int d: {0,1}) {
			c->cd(d+1); 
			h2_tpc_xd[i][d]->Draw("COLZ");
		}
		c->cd(3);
		h2_tpc_xy[i]->Draw("COLZ");
		for(int a: {0,1,2,3}) {
			c->cd(5+a);
			h2_tpc_yd[i][a]->Draw("COLZ");
		}
		c->cd(4);
		PLatex( 0.08, 
			Form("TPC%s", label[i]),
			"Reference made by: ",
			eval_latex_label(0),
			eval_latex_label(1),
			eval_latex_label(2),
			eval_latex_label(3)
		);
	}

	if(e == OptionE::FitDiff) {
		/* Idea is to take central bin positions (xi),
		 * their TH1D* projection has a peak around value (yi) with counts (wi)
		 * and then do a weighted linear fit. */
		double dist_side = 15.0; 
		for(int i=0; i<N; ++i) {
			if(I >= 0 and i != I) continue;

			const TPCParam& p = tpc_params->at(i);
			WARN("\n" EMPH1(TPC%s) " " EBOLD(fitting parameters) "\n", label[i]);
			for(int d: {0,1}) {
				const double b0 = p.x_offset[d];
				const double a0 = p.x_factor[d];

				auto [k,l] = FitProfile(*h2_tpc_xd[i][d], dist_side);
				WARN("DL%d found \'more optimized\' parameter\n", d);
				WARN("\rOffset/slope: %.9f, %.9f\n", l, k);
				printf("\rBefore: (%.6f , %.6f)\n",             b0,                 a0);
				printf("\rNew   : " EBOLD((%.6f , %.6f)) "\n",  b0/(k+1) - l/(k+1), a0/(k+1));
			}
			dist_side = 20.0;
			for(int a: {0,1,2,3}) {
				const double b0 = p.y_offset[a];
				const double a0 = p.y_factor[a];

				auto [k,l] = FitProfile(*h2_tpc_yd[i][a], dist_side);
				WARN("A %d found \'more optimized\' parameter\n", a);
				WARN("\rOffset/slope: %.5f, %.5f\n", l, k);
				printf("\rBefore: (%.6f , %.6f)\n",             b0,                 a0);
				printf("\rNew   : " EBOLD((%.6f , %.6f)) "\n",  b0/(k+1) - l/(k+1), a0/(k+1));
			}
		}
	}
}

double lo_y(TH2D* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetYaxis()->GetXmax() + (1+r)/2 * h2->GetYaxis()->GetXmin(); 
}
double hi_y(TH2D* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetYaxis()->GetXmax() + (1-r)/2 * h2->GetYaxis()->GetXmin(); 
}

std::array<double, 2> FitProfile (
	TH2D* h2, int N_BINS_SIDE
) {
	std::vector<double> xv, yv, wv;
	xv.reserve(2*N_BINS_SIDE + 1);
	yv.reserve(2*N_BINS_SIDE + 1);
	wv.reserve(2*N_BINS_SIDE + 1);
	TAxis* ax = h2->GetXaxis();

	int b0 = ax->FindBin(0.0);
	const int bMin = std::max(             1, b0 - N_BINS_SIDE);
	const int bMax = std::min(ax->GetNbins(), b0 + N_BINS_SIDE);
	for(int bx=bMin; bx <= bMax; ++bx) {
		double x = ax->GetBinCenter(bx);
		auto py = std::unique_ptr<TH1D>( h2->ProjectionY("_py_bx", bx, bx) ); 
		py->SetDirectory(nullptr);
		int by =  py->GetMaximumBin();
		double w = py->GetBinContent(by); 
		double y = py->GetBinCenter(by);
		xv.push_back(x); 
		yv.push_back(y); 
		wv.push_back(w);
	}

	WARN("Doing fit with: ");
	std::cerr << std::endl
		<< "x: " << KBH_BLU << xv << KNRM << std::endl
		<< "y: " << KBH_RED << yv << KNRM << std::endl
		//<< std::endl << "w: " << KBH_MAG << wv << KNRM << std::endl
		;
	//return PolyFit<1>(xv, yv, wv);
	return PolyFit<1>(xv, yv);
}

/* No checks, fuck it. */
std::array<double, 2> FitProfile (
	TH2D* h2, double distance
) {
	TAxis* ax = h2->GetXaxis();
	int ibin_side = ax->FindBin(distance + 0.0000001);
	int b0 = ax->FindBin(0.0);
	
	int bins_to_the_side = std::max(ibin_side - b0, 1);
	return FitProfile(h2, bins_to_the_side);
}
