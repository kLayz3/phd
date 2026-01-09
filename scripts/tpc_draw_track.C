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

std::array<double, 2> MeanXY(const std::array<RNTPCCal::Measurement, 4>& );

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
	constexpr bool _take[][N] = {
		{1, 1, 1, 1}, // 21
		{1, 1, 1, 1}, // 22
		{1, 1, 1, 1}, // 23
		{0, 0, 0, 0}  // 24
	};

	constexpr double WIDTH = 70.0;
	constexpr double ANODE_SEP = 14.0; 
	const Arr2<double, N, 4> z0 = [tpc_param](){
		Arr2<double,N,4> z{};
		for(int i=0; i<N; ++i) {
			z[i][0] = tpc_param->at(i).z0 - WIDTH/2 - ANODE_SEP/2;
			z[i][1] = tpc_param->at(i).z0 - WIDTH/2 + ANODE_SEP/2;
			z[i][2] = tpc_param->at(i).z0 + WIDTH/2 - ANODE_SEP/2;
			z[i][3] = tpc_param->at(i).z0 + WIDTH/2 + ANODE_SEP/2;
		}
		return z;
	}();

	const double zT = (3355 - 440/2); /* target position. */
	

	TH2I* h2_tpc_xy[N];
	for(int i=0; i<N; ++i) 
		h2_tpc_xy[i] = new TH2I(Form("h2_tpc%d_xy", i), Form("TPC%d positions Y vs. X", i), 400, -20, 20, 400, -20, 20); 
	TH2I* h2_xy = new TH2I("h2_xy", "XY-pos at target", 400, -50, 50, 400, -50, 50);
	TH2I* h2_ab = new TH2I("h2_ab", "XY-angle [mrad]", 400, -50, 50, 400, -50, 50);
	TH2I* h2_track_x = new TH2I("h2_track_x", "Track density (X)", 600, 0, 4500, 500, -60, 60); 
	TH2I* h2_track_y = new TH2I("h2_track_y", "Track density (Y)", 600, 0, 4500, 500, -60, 60); 

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	
	bool take_event = true;
	
	std::vector<double> x; x.reserve(4*4);
	std::vector<double> y; y.reserve(4*4);
	std::vector<double> z; x.reserve(4*4);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		x.clear(); y.clear(); z.clear();

		for(int i = 0; i < N; ++i) {
			const auto& tpc = frs->tpc[i];
			if(tpc.hits.size() != 1) continue;

			const std::array<Measurement, 4>& hit = tpc.hits[0];
			for(int j=0; j<4; ++j) {
				if(!std::isnan(hit[j].x) and _take[i][j]) {
					x.push_back( hit[j].x );
					y.push_back( hit[j].y );
					z.push_back( z0[i][j] );
				}
			}
			auto m = MeanXY(hit);
			h2_tpc_xy[i]->Fill( m[0], m[1] );
		}
		if(z.size() < 7) continue;
		auto fx = PolyFit<1>(z, x);	
		auto fy = PolyFit<1>(z, y);	
		
		h2_xy->Fill(fx[1]*zT + fx[0], fy[1]*zT + fy[0]);
		h2_ab->Fill(fx[1]*1000.0, fy[1]*1000);

		FillTrack(h2_track_x, fx);
		FillTrack(h2_track_y, fy);
	}

	TCanvas* cT = new TCanvas("cT", "cT", 2000, 1200);
	cT->Divide(1,2);
	cT->cd(1);
	h2_xy->Draw("COLZ");
	cT->cd(2);
	h2_ab->Draw("COLZ");

	TCanvas* cTr = new TCanvas("cTr", "cTr", 2000, 1200);
	cTr->Divide(1,2);
	cTr->cd(1);
	h2_track_x->Draw("COLZ");
	cTr->cd(2);
	h2_track_y->Draw("COLZ");

	TCanvas* c[N];
	for(int i=0; i<N; ++i) {
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 2000, 1200);
		h2_tpc_xy[i]->Draw("COLZ");
	}
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
