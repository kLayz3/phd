#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooPolynomial.h"

using namespace ROOT;
using namespace ROOT::Experimental;

void tpc_xy_explo(std::string fileName = "", int i = 0, 
	std::array<double, 2> x0 = {0,0},      std::array<double, 2> wx = {6,6},
	std::array<double, 4> y0 = {0,0,0,0}, std::array<double, 4> wy = {6,6,6,6})

{
	//ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH1I* h1_tpc_x[2];
	TH1I* h1_tpc_y[4];

	for(int j=0; j<2; ++j) {
		h1_tpc_x[j] = new TH1I(Form("TPC%d_X%d", i, j), Form("TPC%d x (dline %d)", i, j), 1000, -50, 50);
		h1_tpc_x[j]->GetXaxis()->SetTitle("mm");
		h1_tpc_x[j]->GetYaxis()->SetTitle("count");
		h1_tpc_x[j]->SetFillColor(2+i);
		h1_tpc_x[j]->SetLineColor(kBlack);
	}
	for(int j=0; j<4; ++j) {
		h1_tpc_y[j] = new TH1I(Form("TPC%d_Y%d", i, j), Form("TPC%d y (anode %d)", i, j), 888, -50, 50);
		h1_tpc_y[j]->GetXaxis()->SetTitle("mm");
		h1_tpc_y[j]->GetYaxis()->SetTitle("count");
		h1_tpc_y[j]->SetFillColor(2+i);
		h1_tpc_y[j]->SetLineColor(kBlack);
	}

	std::vector<RooRealVar> xs; xs.reserve(2);
	std::vector<RooDataSet> val_xs; val_xs.reserve(2);
	std::vector<RooRealVar> ys; ys.reserve(2);
	std::vector<RooDataSet> val_ys; val_ys.reserve(2);
	std::vector<RooRealVar> x_mean; x_mean.reserve(2);
	std::vector<RooRealVar> x_sigma; x_sigma.reserve(2);
	std::vector<RooGaussian> x_gauss; x_gauss.reserve(2);

	for(int j=0; j<2; ++j) {
		xs.     emplace_back(Form("x%d",j), Form("x%d",j), x0[j] - 2*wx[j], x0[j] + 2*wx[j]);
		val_xs .emplace_back(Form("val_x%d",j), Form("val_x%d",j), RooArgSet(xs[j])); 
		x_mean .emplace_back(Form("x_mean%d" , j), Form("x_mean%d" , j), x0[j], x0[j] - wx[j], x0[j] + wx[j]);
		x_sigma.emplace_back(Form("x_sigma%d", j), Form("x_sigma%d", j), wx[j]/3, 0.01, wx[j]);
		x_gauss.emplace_back(Form("x_gauss%d", j), Form("x_gauss%d", j), xs[j], x_mean[j], x_sigma[j]);
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		const auto& tpc = frs->tpc[i];
		if(tpc.hits.size() != 1) continue;
		const auto& hit_a4 = tpc.hits[0];
		
		for(int a=0; a<4; ++a) {
			const auto& hit = hit_a4[a];
			int dl = (a >> 1);
			
			auto& fitX = xs[dl];
			auto& fitV = val_xs[dl];
			fitX = hit.x;
			fitV.add(RooArgSet(fitX));

			h1_tpc_x[dl] -> Fill(hit.x);
			h1_tpc_y[a]  -> Fill(hit.y);
		}
	}

/*
	for(int j=0; j<2; ++j) {
		auto& data = val_xs[j];
		auto& gauss = x_gauss[j];
		auto& mean = x_mean[j]; auto& sigma = x_sigma[j];

		RooFitResult* res = gauss.fitTo(data, RooFit::Save(true));

		printf("\n(X)[%d] Mean : %.2f ± %.2f\n", j, mean.getVal(), mean.getError());
		printf("(X)[%d] Sigma: %.2f ± %.2f\n", j, sigma.getVal(), sigma.getError());
		//res->Print("v");
	}
*/
	TCanvas* cX = new TCanvas(Form("TPC%d-X", i), Form("TPC%d-X", i), 1600, 1200);
	cX->Divide(1,2);

	cX->cd(1);
	h1_tpc_x[0]->Draw();

	cX->cd(2);
	h1_tpc_x[1]->Draw();

	TCanvas* cY = new TCanvas(Form("TPC%d-Y", i), Form("TPC%d-Y", i), 1600, 1200);
	cY->Divide(2,2);
	
	cY->cd(1);
	h1_tpc_y[0]->Draw();
	
	cY->cd(2);
	h1_tpc_y[1]->Draw();
	
	cY->cd(3);
	h1_tpc_y[2]->Draw();
	
	cY->cd(4);
	h1_tpc_y[3]->Draw();
}
