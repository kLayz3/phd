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

template<size_t Size>
using Vec = Eigen::Matrix<double, Size, 1>;

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

void tpc_draw_track(std::string fileName = "") {
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
		{1, 1}  // 24
	};
	
	constexpr double WIDTH = 70.0;
	const Arr2<double, N, 2> z0 = [tpc_param](){
		Arr2<double,N,2> z{};
		for(int i=0; i<N; ++i) {
			z[i][0] = tpc_param->at(i).z0 - WIDTH/4;
			z[i][1] = tpc_param->at(i).z0 - WIDTH/4;
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

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	 
	bool take_event = true;
	
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

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		x.clear(); y.clear(); z.clear();
		/* Get referent hit for correlation plot. 
		 * Focus is at roughly ~Z_FOC, extrapolate position from this reference and focus. */
		xref = (frs->tpc[0].hits[0].size() > 0) ? frs->tpc[0].hits[0].at(0).x : NAN;

		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];

			const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;

			for(int d=0; d<2; ++d) {
				const std::vector<Measurement>& hits = tpc_hits[d];
				
				if(hits.size() == 1 and !std::isnan(hits[0].x) and _take[i][d]) {
					x.push_back( hits[0].X() );
					y.push_back( hits[0].Y() );
					z.push_back( z0[i][d] );
				}
			}

			if(tpc_hits[0].size() == 1 and tpc_hits[1].size() == 1) { 
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
		}

		if(z.size() < 3) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	
		
		h2_xy->Fill(fx[1]*zT + fx[0], fy[1]*zT + fy[0]);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		FillTrack(h2_track_x, fx);
		FillTrack(h2_track_y, fy);
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

	TCanvas* c[N];
	for(int i=0; i<N; ++i) {
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 2000, 1200);
		h2_tpc_xy[i]->Draw("COLZ");
	}
}

std::array<double, 2> MeanXY(const std::array<std::vector<RNTPCCal::Measurement>, 2>& tpc_hits) {
	int N = 0;
	std::array<double, 2> p {0,0};
	for(const auto& dl : tpc_hits) {
		if(dl.size() == 1 and !std::isnan(dl[0].x)) {
			p[0] += dl[0].X();
			p[1] += dl[0].Y();
			++N;
		}
	}
	p[0] /= N;
	p[1] /= N;
	return p;
}

double lo_x(TH2I* h2, double r) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetXaxis()->GetXmax() + (1+r)/2 * h2->GetXaxis()->GetXmin(); 
}
double hi_x(TH2I* h2, double r) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetXaxis()->GetXmax() + (1-r)/2 * h2->GetXaxis()->GetXmin(); 
}

double lo_y(TH2I* h2, double r) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetYaxis()->GetXmax() + (1+r)/2 * h2->GetYaxis()->GetXmin(); 
}
double hi_y(TH2I* h2, double r) { 
	if(r > 1 || r <= 0)
		throw std::runtime_error("second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetYaxis()->GetXmax() + (1-r)/2 * h2->GetYaxis()->GetXmin(); 
}
