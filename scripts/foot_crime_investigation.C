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

static constexpr double WIDTH = 440;
static constexpr double z0 = 3355;
static constexpr double zT = z0 - WIDTH/2;
static constexpr double z[] = {
	zT + 89,
	zT + 89 + 90,
	zT + 89 + 90 + 136,
	zT + 89 + 90 + 136 + 100
};
static constexpr double zTPC = 3982;

#if 1 
#	define EXTRAPOLATE_TPC24
#else
#	define EXTRAPOLATE_UPSTREAM

constexpr int N = 3;
constexpr bool _take[N][2] = {
	{1, 1}, // 21
	{1, 1}, // 22
	{0, 0}, // 23
};
constexpr double WIDTH_TPC = 70.0;

#endif

void foot_crime_investigation(std::string fileName = "", uint32_t pos = 0) {
	if(pos > 3) throw std::runtime_error("Second argument `pos` can be only {0,1,2,3}.");
	ROOT::EnableImplicitMT();
	
	using Measurement = RNTPCCal::Measurement;

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	std::array <
		std::shared_ptr<RNFOOTCal>, 8
	> foot {};
	for(int i=0; i<8; ++i) {
		foot[i] = model->MakeField<RNFOOTCal>(Form("FOOT%d", i));
	}
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH2I *h2_foot_x[8];
	TH2I *h2_foot_y[8];
	for(int i=0; i<8; ++i) {
		TH2I* hx = new TH2I(Form("FOOT%d_x_corr", i), Form("Possible: TPC vs. FOOT%d X", i), 320,0,640, 80, -20, 20);
		TH2I* hy = new TH2I(Form("FOOT%d_y_corr", i), Form("Possible: TPC vs. FOOT%d Y", i), 320,0,640, 80, -20, 20);
		hx->GetXaxis()->SetTitle("FOOT strip #");
		hx->GetYaxis()->SetTitle("TPC X-extrapolated [mm]");
		hy->GetXaxis()->SetTitle("FOOT strip #");
		hy->GetYaxis()->SetTitle("TPC Y-extrapolated [mm]");
		h2_foot_x[i] = hx; 
		h2_foot_y[i] = hy; 
	}

#ifdef EXTRAPOLATE_TPC24
	TH2I* h2_tpc = new TH2I("h2_tpc", "TPC y vs. x", 160, -80, 80, 160, -80, 80);
	const double k = ( z[ static_cast<int>(pos) ] - zT ) / (zTPC - zT);

#else	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_param = f->Get < 
			std::remove_reference_t<decltype(*tpc_param)>
		> ("FRS_tpc_parameters");
		if(!tpc_param)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
	}
	const std::array<double, N> z0 = {
		tpc_param->at(0).z0,
		tpc_param->at(1).z0,
		tpc_param->at(2).z0
	};

	std::vector<double> x; x.reserve(N);
	std::vector<double> y; y.reserve(N);
	std::vector<double> z; x.reserve(N);
#endif

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		double x_extrapolated, y_extrapolated;

#ifdef EXTRAPOLATE_TPC24
		const auto& tpc24 = frs->tpc.at(3);
		double x24 = tpc24.X0();
		double y24 = tpc24.Y0();

		h2_tpc->Fill(x24, y24);

		x_extrapolated = k * x24; 
		y_extrapolated = k * y24; 
#else
		x.clear(); y.clear(); z.clear();
		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];

			double x0 = tpc.X0();
			double y0 = tpc.Y0();
			
			if(!std::isnan(x0) and !std::isnan(y0)) {
				x.push_back(x0);
				y.push_back(y0);
				z.push_back(z0[i]);
			}		
		}
		
		if(z.size() < 2) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	
		
		x_extrapolated = fx[0] + fx[1]*z[pos]; 
		y_extrapolated = fy[0] + fy[1]*z[pos]; 

#endif
		for(int i=0; i<8; ++i) {
			std::vector<double> cl_pos = foot[i]->X();
			for(const auto p : cl_pos) {
				h2_foot_x[i]->Fill(p, x_extrapolated);
				h2_foot_y[i]->Fill(p, y_extrapolated);
			}
		}
	}

	TCanvas* c[8];
	for(int i=0; i<8; ++i) {
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 1600, 1200);
		c[i]->Divide(1,2);
		c[i]->cd(1);
		gPad->SetGrid();
		h2_foot_x[i]->Draw("COLZ");
		c[i]->cd(2);
		gPad->SetGrid();
		h2_foot_y[i]->Draw("COLZ");
	}

#ifdef EXTRAPOLATE_TPC24
	TCanvas* c_tpc = new TCanvas("c_tpc", "c_tpc", 1600, 1200);
	h2_tpc->Draw("COLZ");

	printf("Positions at FOOT taken from extrapolating TPC24 measurement!\n");
#else
	printf("Positions at FOOT taken from extrapolating TPC21/22/23 track!\n");
#endif

}
