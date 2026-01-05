#include "libs.hh"
#include "Eigen/Dense"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;

template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

template<size_t Size>
using Vec = Eigen::Matrix<double, Size, 1>;

int main(int argc, char* argv[]) {
	//ROOT::EnableImplicitMT();
	TApplication* app = new TApplication("myApp", 0, 0);
	
	if(argc < 2) {
		perror("Supply file name");
		exit(1);
	}

	const char* fileName = argv[1];

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName, "READ");
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

	const double* z = &z0[0][0]; // Used to iterate.
	const double zT = (3355 - 440.0/2);
	
	TH1I* h1_residual[N];
	TH1I* h1_x_target, *h1_a;

	for(int i=0; i<N; ++i) {
		h1_residual[i] = new TH1I(Form("h1_resi_%d", i), Form("TPC-2%d [%d] residual", i/2 + 1, i%2 + 1),
			200, -100, 100);
	}
	h1_x_target = new TH1I("h1_x_target", "X-pos at target", 200, -50, 50);
	h1_a = new TH1I("h1_a", "X-angle [mrad]", 200, -50, 50);

	const Eigen::Matrix<double,N,2> A = [&] {
		Eigen::Matrix<double,N,2> A;
		for(int i=0; i<N; ++i) {
			if(take[i]) {
				A(i,0) = z[i];
				A(i,1) = 1.0;
			} else {
				A(i,0) = 0;
				A(i,1) = 0;
			}
		} 
		return A;
	} ();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	using Measurement = RNTPCCal::Measurement;
	
	bool take_event = true;
	
	Vec<N> x;
	
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		take_event = true;
		for(int itpc = 0; itpc<4; ++itpc) {
			for(int idl=0; idl<2; ++idl ) {
				int i = itpc*2 + idl;
				if(!take[i]) continue;
			
				const auto& tpc = frs->tpc[itpc];
				if(tpc.hits.size() != 0) take_event = false;
				const Measurement& m = tpc.hits[0][idl * 2];
				
				if(std::isnan(m.x)) take_event = false; 
				x(i) = m.x; 
			}
		}
		
		if(!take_event) continue;

		Vec<2> fit = A.colPivHouseholderQr().solve(x);
		
		double a   = fit(0);
		double off = fit(1);
			
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

	printf("End of main."); 
	app->Run();
	return 0;
}
