#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../includes/PolyFitter.hxx"
#include "../includes/Tracking.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

std::array<double, 2> MeanXY(const std::array<std::vector<RNTPCCal::Measurement>, 2>& tpc_hits);

static const char* labels[] = {
	"21", "22", "23", "24", "31", "41", "42"
};
constexpr double zT = 3355 - 440/2;

double lo_x(TH2I*, double );
double hi_x(TH2I*, double );
double lo_y(TH2I*, double );
double hi_y(TH2I*, double );

//#define PRINT_ME
using Suspect = std::pair<int, int>;
bool has_suspect(Suspect& );

constexpr std::pair<int,int> combo[] {
	{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}
};

void tpc_draw_track(std::string fileName = "", std::pair<int, int> sus = {-1,-1}, uint64_t max_events = -1 ) {
	if(has_suspect(sus)) {
		if(sus.first < 0 or sus.first > 3) throw std::runtime_error("Suspect tpc id must be {0,1,2,3}.");
		if(sus.second < 0 or sus.second > 1) throw std::runtime_error("Suspect delay line must be {0,1}.");
	}

	ROOT::EnableImplicitMT();
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
	}

	constexpr int N = 4;
	constexpr bool _take[N][2] = {
		{1, 1}, // 21
		{1, 1}, // 22
		{1, 1}, // 23
		{0, 0}  // 24
	};
	
	constexpr double WIDTH = 60.0;
	const Arr2<double, N, 2> zDL = [tpc_param](){
		Arr2<double,N,2> z{};
		for(int i=0; i<N; ++i) {
			z[i][0] = tpc_param->at(i).z0 - WIDTH/4;
			z[i][1] = tpc_param->at(i).z0 + WIDTH/4;
		}
		return z;
	}();
	const std::array<double, N> zTPC = {
		tpc_param->at(0).z0,	
		tpc_param->at(1).z0,	
		tpc_param->at(2).z0,	
		tpc_param->at(3).z0
	};
	printf("TPC positions: \n");
	for(int i=0; i<N; ++i) printf("TPC%s: %.1f mm\n", labels[i], zTPC[i]);

	TH2I* h2_tpc_xy[N];
	for(int i=0; i<N; ++i) {
		h2_tpc_xy[i] = new TH2I(Form("h2_tpc%d_xy", i), Form("TPC%s positions Y vs. X", labels[i]), 400, -20, 20, 400, -20, 20); 
		h2_tpc_xy[i]->GetXaxis()->SetTitle("x [mm]");
		h2_tpc_xy[i]->GetYaxis()->SetTitle("y [mm]");
	}

	TH2I* h2_xy = new TH2I("h2_xy", "XY-position at target", 400, -50, 50, 400, -50, 50);
	h2_xy->GetXaxis()->SetTitle("x [mm]"); h2_xy->GetYaxis()->SetTitle("x [mm]");
	TH2I* h2_ab = new TH2I("h2_ab", "XY-angle at target", 400, -50, 50, 400, -50, 50);
	h2_ab->GetXaxis()->SetTitle("a[x] [mrad]"); h2_ab->GetYaxis()->SetTitle("b[y] [mrad]");

	TH2I* h2_track_x = new TH2I("h2_track_x", "Track density S2 (X)", 600, 0, 4500, 500, -60, 60); 
	TH2I* h2_track_y = new TH2I("h2_track_y", "Track density S2 (Y)", 600, 0, 4500, 500, -60, 60); 
	h2_track_x->GetXaxis()->SetTitle("z [mm]"); h2_track_x->GetYaxis()->SetTitle("x [mm]");
	h2_track_y->GetXaxis()->SetTitle("z [mm]"); h2_track_y->GetYaxis()->SetTitle("y [mm]");
	
	constexpr int N_BINS_X  = 100;
	constexpr int PROJ_X_LO = -50;
	constexpr int PROJ_X_HI =  50;
	constexpr int N_BINS_Y  = 100;
	constexpr int PROJ_Y_LO = -20;
	constexpr int PROJ_Y_HI =  20;

	int N0, N1;
	
	TH2I* h2_tpc_xproj[4][3];
	TH2I* h2_tpc_yproj[4][3];
	TH2I* h2_tpc_sus_xproj[2][2];
	TH2I* h2_tpc_sus_yproj[2][2];
	
	if(has_suspect(sus)) {
		printf("\nInvestigating suspicious TPC%s, delay line: %d\n", labels[sus.first], sus.second);
		for(int i=0; i<4; ++i) {
			if(i == sus.first) continue;
			printf("%d => TPC%s\n", i, labels[i]);
		}
		std::cout << "TPC" << labels[sus.first] << " delay line " << sus.second << ";  enter two TPC indices to take as ref: ";
		std::cin >> N0;
		if(N0 != sus.first and N0 >= 0 and N0 < 4)
			std::cout << "Fine, seen: " << labels[N0] << " as initial reference.\n";
		else 
			throw std::runtime_error(Form("Invalid TPC ID, must be {0,1,2,3} and different than: %d\n", sus.first));
		std::cout << "Enter second TPC index: ";
		std::cin >> N1;
		if(N1 != sus.first and N1 >= 0 and N1 < 4 and N1 != N0)
			std::cout << "Fine, seen: " << labels[N1] << " as second reference. Proceeding with analysis ...\n";
		else
			throw std::runtime_error(Form("Invalid TPC ID, must be {0,1,2,3} and different than: %d, %d\n", sus.first, N0));
	
		for(int d0: {0,1}) {
			for(int d1: {0,1}) {
				h2_tpc_sus_xproj[d0][d1] = new TH2I(Form("h2_tpc_sus_xproj%d-%d", d0, d1), 
					Form("TPC%s(%d) && TPC%s(%d) X-projection diff on TPC%s(%d)", labels[N0],d0,labels[N1],d1,labels[sus.first],sus.second),
						N_BINS_X, PROJ_X_LO, PROJ_X_HI, 60, -10, 10);
				h2_tpc_sus_xproj[d0][d1]->GetXaxis()->SetTitle(Form("TPC%s(%d) X-measurement [mm]", labels[sus.first], sus.second)); 
				h2_tpc_sus_xproj[d0][d1]->GetYaxis()->SetTitle(Form("TPC-X (%s:%d && %s:%d) extrap. - TPC%s(%d) X [mm]", 
					labels[N0],d0, labels[N1],d1, labels[sus.first],sus.second)); 
				
				h2_tpc_sus_yproj[d0][d1] = new TH2I(Form("h2_tpc_sus_yproj%d-%d", d0, d1), 
					Form("TPC%s(%d) && TPC%s(%d) Y-projection diff on TPC%s(%d)", labels[N0],d0,labels[N1],d1,labels[sus.first],sus.second),
						N_BINS_Y, PROJ_Y_LO, PROJ_Y_HI, 60, -10, 10);
				h2_tpc_sus_yproj[d0][d1]->GetXaxis()->SetTitle(Form("TPC%s(%d) Y-measurement [mm]", labels[sus.first], sus.second)); 
				h2_tpc_sus_yproj[d0][d1]->GetYaxis()->SetTitle(Form("TPC-Y (%s:%d && %s:%d) extrap. - TPC%s(%d) Y [mm]", 
					labels[N0],d0, labels[N1],d1, labels[sus.first],sus.second)); 
			}
		}
	} 
	else {
		for(int itpc: {0,1,2,3}) {
			/* 3 possible combinations that could extrapolate into this `itpc` index of TPC. */
			int i = 0;
			for(auto [i0,i1]: combo) {
				char name[13] = {'\0'};
				if(i0 == itpc || i1 == itpc) continue;
				snprintf(name, 13, "%s&%s => %s", labels[i0], labels[i1], labels[itpc]);


				h2_tpc_yproj[itpc][i] = new TH2I(Form("h2_tpc_yproj%d-%d",  itpc, i), Form("TPC: %s Y-projection diff", name),
						N_BINS_X, PROJ_X_LO, PROJ_X_HI, 60, -10, 10);
				h2_tpc_xproj[itpc][i] = new TH2I(Form("h2_tpc_xproj%d-%d",  itpc, i), Form("TPC: %s X-projection diff", name),
						N_BINS_Y, PROJ_Y_LO, PROJ_Y_HI, 60, -10, 10);
				h2_tpc_xproj[itpc][i]->GetXaxis()->SetTitle(Form("TPC%s X-measurement [mm]", labels[itpc])); 
				h2_tpc_xproj[itpc][i]->GetYaxis()->SetTitle(Form("TPC X (%s&%s) - TPC%s X-measurement [mm]", labels[i0], labels[i1], labels[itpc])); 

				h2_tpc_yproj[itpc][i]->GetXaxis()->SetTitle(Form("TPC%s Y-measurement [mm]", labels[itpc])); 
				h2_tpc_yproj[itpc][i]->GetYaxis()->SetTitle(Form("TPC Y (%s&%s) - TPC%s Y-measurement [mm]", labels[i0], labels[i1], labels[itpc])); 
				++i;
			}
		}
	}

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	 
	std::vector<double> x; x.reserve(4*2);
	std::vector<double> y; y.reserve(4*2);
	std::vector<double> z; x.reserve(4*2);
	
#ifdef PRINT_ME
	uint32_t printme[4] = {0}; constexpr double MAX_PRINT = 10;
	constexpr double prob = 0.1;
#endif
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

			auto m = MeanXY(tpc_hits);
			h2_tpc_xy[i]->Fill( m[0], m[1] );
		}

		if(z.size() < 3) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	

		h2_xy->Fill(fx[1]*zT + fx[0], fy[1]*zT + fy[0]);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		FillTrack(h2_track_x, fx);
		FillTrack(h2_track_y, fy);

		if(has_suspect(sus)) {
			auto [NR, ndl] = sus;
			
			if(frs->tpc[NR].hits[ndl].size() == 0) continue;;
			double xR = frs->tpc[NR].hits[ndl][0].X();
			double yR = frs->tpc[NR].hits[ndl][0].Y();
			//double zR = zDL[NR][ndl];
			double zR = zTPC[NR];
			if(std::isnan(xR) || std::isnan(yR)) continue;
				
			for(int d0: {0,1}) {
				if(frs->tpc[N0].hits[d0].size() == 0) continue;
				double x0 = frs->tpc[N0].hits[d0][0].X();
				double y0 = frs->tpc[N0].hits[d0][0].Y();
				//double z0 = zDL[N0][d0];
				double z0 = zTPC[N0];
				
				for(int d1: {0,1}) {
					if(frs->tpc[N1].hits[d1].size() == 0) continue;
					double x1 = frs->tpc[N1].hits[d1][0].X();
					double y1 = frs->tpc[N1].hits[d1][0].Y();
					//double z1 = zDL[N1][d1];
					double z1 = zTPC[N1];

					double kx = (x1 - x0) / (z1 - z0);
					double ky = (y1 - y0) / (z1 - z0);
					h2_tpc_sus_xproj[d0][d1] -> Fill (
						xR,
						(kx * (zR - z0) + x0) - xR
					);
					h2_tpc_sus_yproj[d0][d1] -> Fill (
						yR,
						(ky * (zR - z0) + y0) - yR
					);
				}
			}
		}
		else { /* Extrapolations per TPC, only if `suspect` is not really given. */
			for(int iref = 0; iref < 4; ++iref) {
				int i=0;
				for(auto [i0,i1]: combo) {
					if(i0 == iref || i1 == iref) continue;

					double xR = frs->tpc[iref].X0();
					double yR = frs->tpc[iref].Y0();
					double x0 = frs->tpc[i0].X0();
					double y0 = frs->tpc[i0].Y0();
					double x1 = frs->tpc[i1].X0();
					double y1 = frs->tpc[i1].Y0();
					if(std::isnan(xR) || std::isnan(yR)) continue;
					if(std::isnan(x0) || std::isnan(y0)) continue;
					if(std::isnan(x1) || std::isnan(y1)) continue;

					double zR = zTPC[iref];
					double z0 = zTPC[i0];
					double z1 = zTPC[i1];

					double kx = (x1 - x0) / (z1 - z0);
					double ky = (y1 - y0) / (z1 - z0);
					h2_tpc_xproj[iref][i] -> Fill (
						xR,
						(kx * (zR - z0) + x0) - xR
					);
					h2_tpc_yproj[iref][i] -> Fill (
						yR,
						(ky * (zR - z0) + y0) - yR
					);
					++i;
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
		TLine* line = new TLine( zTPC[i], lo_y(h2_track_x, r), zTPC[i], hi_y(h2_track_x, r) );
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	TLine* line = new TLine( zT, lo_y(h2_track_x, r), zT, hi_y(h2_track_x, r) );
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");

	cTr->cd(2);
	h2_track_y->Draw("COLZ");
	for(int i=0; i<4; ++i) {
		TLine* line = new TLine( zTPC[i], lo_y(h2_track_y, r), zTPC[i], hi_y(h2_track_y, r) );
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		line->Draw("SAME");
	}
	line = new TLine( zT, lo_y(h2_track_y, r), zT, hi_y(h2_track_y, r) );
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	line->Draw("SAME");

	/* ========================= */

#if 0
	TCanvas* c[N];
	for(int i=0; i<N; ++i) {
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 2000, 1200);
		h2_tpc_xy[i]->Draw("COLZ");
	}
#endif
	
	if(has_suspect(sus)) {
		auto [NR, ndl] = sus;
		TCanvas* c = new TCanvas(Form("cP%s", labels[NR]), Form("Corr: TPC%s-%d", labels[NR], ndl), 1900, 1400);
		c->Divide(4,2); int cid = 1;
		for(int d0: {0,1}) {
			for(int d1: {0,1}) {
				c->cd( cid++ );
				gPad->SetGrid();
				h2_tpc_sus_xproj[d0][d1]->Draw("COLZ");
			}
		} 
		for(int d0: {0,1}) {
			for(int d1: {0,1}) {
				c->cd( cid++ );
				gPad->SetGrid();
				h2_tpc_sus_yproj[d0][d1]->Draw("COLZ");
			}
		}
	}
	else {
		for(int itpc: {0,1,2,3}) {
			TCanvas* c = new TCanvas(Form("cP%s", labels[itpc]), Form("cP%s", labels[itpc]), 1800, 1400);
			c->Divide(3,2);
			
			TH2I** h2_tpc_x = h2_tpc_xproj[itpc]; // [3];
			TH2I** h2_tpc_y = h2_tpc_yproj[itpc]; // [3];
			
			for(int i: {0,1,2}) {
				c->cd(i+1);
				gPad->SetGrid();
				h2_tpc_x[i]->Draw("COLZ");

				c->cd(i+4);
				gPad->SetGrid();
				h2_tpc_y[i]->Draw("COLZ");
			}
		}
	}
}

std::array<double, 2> MeanXY(const std::array<std::vector<RNTPCCal::Measurement>, 2>& tpc_hits) {
	int Nx = 0; int Ny = 0;
	std::array<double, 2> p {0,0};
	for(const auto& dl : tpc_hits) {
		if(dl.size() == 1 and !std::isnan(dl[0].x)) {
			p[0] += dl[0].x;
			++Nx;
			if(!std::isnan(dl[0].y[0])) {
				p[1] += dl[0].y[0]; ++Ny;
			}
			if(!std::isnan(dl[0].y[1])) {
				p[1] += dl[0].y[1]; ++Ny;
			}
		}
	}
	p[0] /= Nx;
	p[1] /= Ny;
	return p;
}

bool has_suspect(Suspect& s) {
	return s.first != -1 and s.second != -1;
}

double lo_x(TH2I* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetXaxis()->GetXmax() + (1+r)/2 * h2->GetXaxis()->GetXmin(); 
}
double hi_x(TH2I* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetXaxis()->GetXmax() + (1-r)/2 * h2->GetXaxis()->GetXmin(); 
}

double lo_y(TH2I* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetYaxis()->GetXmax() + (1+r)/2 * h2->GetYaxis()->GetXmin(); 
}
double hi_y(TH2I* h2, double r=0) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetYaxis()->GetXmax() + (1-r)/2 * h2->GetYaxis()->GetXmin(); 
}
