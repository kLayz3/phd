#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../includes/PolyFitter.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

template<size_t Size>
using Vec = Eigen::Matrix<double, Size, 1>;

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
	constexpr bool _take[][2] = {
		{1, 1}, // 21
		{1, 1}, // 22
		{1, 1}, // 23
		{1, 1}  // 24
	};
	constexpr int N = (int)(sizeof _take / sizeof *_take[0]);
	
	printf("N (constexpr) = %d\n", N);
	const bool* take = &_take[0][0];

	constexpr double WIDTH = 70.0;
	const Arr2<double, 4, 2> z0 = [tpc_param] {
		Arr2<double,4,2> z{};

		for(int i=0; i<4; ++i) {
			z[i][0] = tpc_param->at(i).z0 - WIDTH / 4;
			z[i][1] = tpc_param->at(i).z0 + WIDTH / 2;
		}
		return z;
	}();

	const double zT = (3355 - 440/2);
	
	TH1I* h1_residual[N];
	TH1I* h1_x_target, *h1_a;

	for(int i=0; i<N; ++i) {
		h1_residual[i] = new TH1I(Form("h1_resi_%d", i), Form("TPC-2%d [%d] residual", i/2 + 1, i%2 + 1),
			120, -30, 30);
	}
	h1_x_target = new TH1I("h1_x_target", "X-pos at target", 200, -50, 50);
	h1_a = new TH1I("h1_a", "X-angle [mrad]", 200, -50, 50);
	
	std::vector<double> z;
	for(int i=0; i<4; ++i) 
		for(int j=0; j<2; ++j)
			if(_take[i][j]) z.push_back(z0[i][j]);

	StaticPolyFitter<N, 1> fitter (z);
		
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	
	bool take_event = true;
	
	std::vector<double> x; x.reserve(N);
	
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		take_event = true;
		for(int itpc = 0; itpc<4; ++itpc) {
			for(int idl=0; idl<2; ++idl ) {
				int i = itpc*2 + idl;
				if(_take[itpc][idl]) continue;
			
				const auto& tpc = frs->tpc[itpc];
				if(tpc.hits.size() != 0) take_event = false;
				const Measurement& m = tpc.hits[0][idl * 2];
				
				if(std::isnan(m.x)) take_event = false; 
				x.push_back(m.x);
			}
		}
		
		if(!take_event) continue;

		assert((z.size() == x.size(), "Sizes of x and z not matching"));
		auto [off, a] = fitter(x);	
			
		h1_a->Fill(a*1000);
		h1_x_target->Fill(a*zT + off);
		
		for(int itpc = 0; itpc<4; ++itpc) {
			for(int idl=0; idl<2; ++idl ) {
				int i = itpc*2 + idl;
				if(!take[i]) continue;
				const Measurement& m = frs->tpc[itpc].hits[0][idl * 2];
				h1_residual[i]->Fill(
					m.x - (a*z[i] + off)
				);
			}
		}
	}

	TCanvas* cres = new TCanvas("cres", "cres", 2000, 1400);
	cres->Divide(2,4);
	for(int i=0; i<N; ++i) {
		cres->cd(i+1);
		h1_residual[i]->Draw();
	}

	TCanvas* cT = new TCanvas("cT", "cT", 2000, 1200);
	cT->Divide(1,2);
	cT->cd(1);
	h1_x_target->Draw();
	cT->cd(2);
	h1_a->Draw();
}
