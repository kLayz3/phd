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

void tpc_draw_track(std::string fileName = "", uint64_t max_events = -1, int do_tpc_p = -1) {
	if(do_tpc_p > 3) throw std::runtime_error("Third argument `pos` can be only {0,1,2,3} or <0 to be defaulted.");

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
		{0, 0}, // 23
		{0, 0}  // 24
	};
	
	constexpr double WIDTH = 70.0;
	const Arr2<double, N, 2> z0 = [tpc_param](){
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
	
	TH2I* h2_tpc_xcorr[4][2];
	TH1I* h1_tpc_xres[4][2];
	TH1I* h1_tpc_yres[4][4];
	TH2I* h2_tpc_xproj[3];
	TH2I* h2_tpc_yproj[3];
	for(int i=0; i<4; ++i) {
		for(int j=0; j<2; ++j) {
			h2_tpc_xcorr[i][j] = new TH2I(Form("h2_tpc%d_xcorr_%d", i,j), Form("h2_tpc%d_xcorr_%d", i,j),
				100, -50, 50, 100, -50, 50);
			h2_tpc_xcorr[i][j]->GetYaxis()->SetTitle("TPC21 dl(0) [mm]");
			h2_tpc_xcorr[i][j]->GetXaxis()->SetTitle(Form("TPC%s dl(%d) [mm]", labels[i], j));
			h1_tpc_xres[i][j] = new TH1I(Form("h1_tpc%d_xres_%d", i,j), Form("TPC%s X-fit residue. Delay-line %d", labels[i], j), 
				100, -20, 20);
			h1_tpc_xres[i][j]->GetXaxis()->SetTitle("Fit residue [mm]");
		}
		for(int a=0; a<4; ++a) {
			h1_tpc_yres[i][a] = new TH1I(Form("h1_tpc%d_yres_%d", i,a), Form("TPC%s Y-fit residue. Anode %d", labels[i], a), 
				100, -20, 20);
			h1_tpc_yres[i][a]->GetXaxis()->SetTitle("Fit residue [mm]");
		}
	}
	for(int combo_id : {0,1,2}) {
		char name[12] = {'\0'};
		if(combo_id == 0) {
			strcpy(name, "21&22 -> 23");
		} else if(combo_id == 1) {
			strcpy(name, "21&23 -> 22");
		} else if(combo_id == 2) {
			strcpy(name, "22&23 -> 21");
		}

		h2_tpc_yproj[ combo_id ] = new TH2I(Form("h2_tpc_yproj%d", combo_id), Form("TPC%s Y-projection diff", name),
			200, -50, 50, 200, -50, 50);
		h2_tpc_xproj[ combo_id ] = new TH2I(Form("h2_tpc_xproj%d", combo_id), Form("TPC%s X-projection diff", name),
			200, -50, 50, 200, -50, 50);
		char* p0 = strtok(name, " -> ");
		char* p1 = strtok(NULL, " -> ");
		h2_tpc_xproj[ combo_id ]->GetXaxis()->SetTitle(Form("TPC%s X-projection [mm]", p0)); 
		h2_tpc_xproj[ combo_id ]->GetYaxis()->SetTitle(Form("TPC%s X [mm]", p1)); 

		h2_tpc_yproj[ combo_id ]->GetXaxis()->SetTitle(Form("TPC%s Y-projection [mm]", p0)); 
		h2_tpc_yproj[ combo_id ]->GetYaxis()->SetTitle(Form("TPC%s Y [mm]", p1)); 
	}

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	 
	std::vector<double> x; x.reserve(4*2);
	std::vector<double> y; y.reserve(4*2);
	std::vector<double> z; x.reserve(4*2);
	
	constexpr double Z_FOC = 2250;
	const double zref = zTPC[0];
	const double λ = 1.0 / (Z_FOC - zref);
	std::array<double, 4> c_x_corr = {
		λ * (Z_FOC - zTPC[0]), 
		λ * (Z_FOC - zTPC[1]), 
		λ * (Z_FOC - zTPC[2]), 
		λ * (Z_FOC - zTPC[3])
	};
	double xref;
	
#ifdef PRINT_ME
	uint32_t printme[4] = {0}; constexpr double MAX_PRINT = 10;
	constexpr double prob = 0.1;
#endif
	uint64_t maxEntries = std::min( max_events, static_cast<uint64_t>(ntuple->GetNEntries()) );
	
	int total_hit_mask = 0;

	for(uint64_t entryId = 0; entryId < maxEntries; ++entryId) {
		ntuple->LoadEntry(entryId);
		total_hit_mask = 0;
		x.clear(); y.clear(); z.clear();
		/* Get referent hit for correlation plot. 
		 * Focus is at roughly ~Z_FOC, extrapolate position from this reference and focus. */
		xref = (frs->tpc[0].hits[0].size() > 0) ? frs->tpc[0].hits[0].at(0).x : NAN;

		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];

			const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;

			for(int d=0; d<2; ++d) {
				const std::vector<Measurement>& hits = tpc_hits[d];
				
				if(hits.size() == 1 and !std::isnan(hits[0].x) and _take[i][d]
					and !std::isnan(hits[0].y[0]) and !std::isnan(hits[0].y[1])) {
					x.push_back( hits[0].X() );
					y.push_back( hits[0].Y() );
					z.push_back( z0[i][d] );
				}
			}

			if(tpc_hits[0].size() == 1 and tpc_hits[1].size() == 1) { 
				h2_tpc_xcorr[i][0]->Fill( tpc_hits[0][0].x, c_x_corr[i] * xref);
				h2_tpc_xcorr[i][1]->Fill( tpc_hits[1][0].x, c_x_corr[i] * xref);
#ifdef PRINT_ME
				if(printme[i] < MAX_PRINT and ((double)rand()/RAND_MAX < prob)) {
					printf("TPC%s: (%.2f, %.2f, %hhu)\t(%.2f, %2.f, %hhu)\n",
						labels[i], tpc_hits[0][0].x, tpc_hits[0][0].y, tpc_hits[0][0].mask, 
						tpc_hits[1][0].x, tpc_hits[1][0].y, tpc_hits[1][0].mask);
					++printme[i];
				}
#endif
			}

			auto m = MeanXY(tpc_hits);
			h2_tpc_xy[i]->Fill( m[0], m[1] );

			if(i < 3 and !std::isnan(m[0]) and !std::isnan(m[1])) 
				total_hit_mask |= (1 << i);
		}

		if(z.size() < 3) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	

		h2_xy->Fill(fx[1]*zT + fx[0], fy[1]*zT + fy[0]);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		FillTrack(h2_track_x, fx);
		FillTrack(h2_track_y, fy);

		/* Residuals: Only fill them if the point went into the Fit calculation. */
		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];
			const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;

			for(int d : {0,1}) {
				if(!_take[i][d]) continue;

				const std::vector<Measurement>& hits = tpc_hits[d];
				if(hits.size() == 1) {
					h1_tpc_xres[i][d] -> Fill( 
						hits[0].X() - 
						(fx[0] + fx[1] * z0[i][d])
					);

					for(int a : {0,1}) {
						h1_tpc_yres[i][2*d + a] -> Fill( 
							hits[0].y[a] - 
							(fy[0] + fy[1] * z0[i][d])
						);
					}
				}
			}
		}

		/* Extrapolations. Only perform if the `total_hit_mask` is valid. */
		if(total_hit_mask == 0b111) {
			double x1 = frs->tpc[0].hits[0][0].X();
			double y1 = frs->tpc[0].hits[0][0].Y();
			double z1 = zTPC[0];

			double x2 = frs->tpc[1].hits[0][0].X();
			double y2 = frs->tpc[1].hits[0][0].Y();
			double z2 = zTPC[1];

			double x3 = frs->tpc[2].hits[0][0].X();
			double y3 = frs->tpc[2].hits[0][0].Y();
			double z3 = zTPC[2];
			
			double kx, ky;
			
			/* 21&22 -> 23 */
			kx = (x2 - x1) / (z2 - z1);
			ky = (y2 - y1) / (z2 - z1);
			h2_tpc_xproj[0]->Fill (
				kx * (z3 - z1) + x1,
				x3
			);
			h2_tpc_yproj[0]->Fill (
				ky * (z3 - z1) + y1,
				y3
			);
			
			/* 21&23 -> 22 */
			kx = (x3 - x1) / (z3 - z1);
			ky = (y3 - y1) / (z3 - z1);
			h2_tpc_xproj[1]->Fill (
				kx * (z2 - z1) + x1,
				x2
			);
			h2_tpc_yproj[1]->Fill (
				ky * (z2 - z1) + y1,
				y2
			);
			
			/* 22&23 -> 21 */
			kx = (x3 - x2) / (z3 - z2);
			ky = (y3 - y2) / (z3 - z2);
			h2_tpc_xproj[2]->Fill (
				kx * (z1 - z2) + x2,
				x1
			);
			h2_tpc_yproj[2]->Fill (
				ky * (z1 - z2) + y2,
				y1
			);
		}
	}
	
	/* Also draw some nominal indicators where each TPC sits. */
	TCanvas* cT = new TCanvas("cT", "cT", 2000, 1200);
	cT->Divide(1,2);
	cT->cd(1);
	h2_xy->Draw("COLZ");

	cT->cd(2);
	h2_ab->Draw("COLZ");

	/* ========================= */

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

	/* ========================= */
	TCanvas* ccorr = new TCanvas("ccorr", "ccorr", 2000, 1400);
	ccorr->Divide(4,2);
	for(int i=0; i<4; ++i) {
		for(int j=0; j<2; ++j) {
			ccorr->cd(2*i + j + 1);
			gPad->SetGrid();
			h2_tpc_xcorr[i][j]->Draw("COLZ");
		}
	}

	TCanvas* cRes[N];
	for(int i=0; i<N; ++i) {
		TCanvas* c = cRes[i] = new TCanvas(Form("cRes-%s", labels[i]), Form("cRes-%s", labels[i]), 1600, 1200);
		c->Divide(3,2);
		for(int d : {0,1}) {
			c->cd(d+1);
			h1_tpc_xres[i][d]->SetFillColor(kGreen);
			h1_tpc_xres[i][d]->Draw();
		}
		for(int a : {0,1,2,3}) {
			c->cd(3+a);
			h1_tpc_yres[i][a]->SetFillColor(kRed);
			h1_tpc_yres[i][a]->Draw();
		}
	}

	TCanvas* cProj = new TCanvas("cProj", "cProj", 1800, 1400);
	cProj->Divide(3,2);
	for(int i: {0,1,2}) {
		cProj->cd(2*i + 1);
		gPad->SetGrid();
		h2_tpc_xproj[i]->Draw("COLZ");

		cProj->cd(2*i + 2);
		gPad->SetGrid();
		h2_tpc_yproj[i]->Draw("COLZ");
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
