#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

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

std::array<double, 2> MeanXY(const std::array<RNTPCCal::Measurement, 4>& );

void foot_crime_investigation(std::string fileName = "", uint32_t pos = 0) {
	if(pos > 4) throw std::runtime_error("Second argument `pos` can be only {0,1,2,3}.");
	ROOT::EnableImplicitMT();
	
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
		h2_foot_x[i] = new TH2I(Form("FOOT%d_x_corr", i), Form("FOOT%d_x_corr", i), 640,0,640, 200, -40, 40);
		h2_foot_y[i] = new TH2I(Form("FOOT%d_y_corr", i), Form("FOOT%d_y_corr", i), 640,0,640, 200, -40, 40);
	}
	TH2I* h2_tpc = new TH2I("h2_tpc", "TPC y vs. x", 200, -40, 40, 200, -40, 40);

	double k = ( z[ static_cast<int>(pos) ] - zT ) / (zTPC - zT);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		auto& tpc24 = frs->tpc.at(3);
		if(tpc24.hits.size() != 1) continue;

		const auto& m = tpc24.hits[0];

		auto p_tpc = MeanXY(m);
		h2_tpc->Fill( p_tpc[0], p_tpc[1] );

		double x_extrapolated = k * p_tpc[0]; 
		double y_extrapolated = k * p_tpc[1]; 
		
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
		h2_foot_x[i]->Draw("COLZ");
		c[i]->cd(2);
		h2_foot_y[i]->Draw("COLZ");
	}
	
	TCanvas* c_tpc = new TCanvas("c_tpc", "c_tpc", 1600, 1200);
	h2_tpc->Draw("COLZ");
}

std::array<double, 2> MeanXY(const std::array<RNTPCCal::Measurement, 4>& hit) {
	int N = 0;
	std::array<double, 2> p {};
	for(const auto& m : hit) {
		if(!std::isnan(m.x)) {
			p[0] += m.x;
			p[1] += m.y;
			++N;
		}
	}
	p[0] /= N;
	p[1] /= N;
	return p;
}
