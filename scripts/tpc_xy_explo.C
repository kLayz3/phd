#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

void tpc_xy_explo(std::string fileName = "") {
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	//constexpr auto N_TPC = RNFRSCal::N_VALID_TPC;
	constexpr auto N_TPC = 3;

	TH1I* h1_tpc_x[N_TPC][2];
	TH1I* h1_tpc_y[N_TPC][4];

	for(int i=0; i<N_TPC; ++i) {
		for(int j=0; j<2; ++j) {
			h1_tpc_x[i][j] = new TH1I(Form("TPC%d_X%d", i, j), Form("TPC%d x (dline %d)", i, j), 180, -30, 30);
			h1_tpc_x[i][j]->GetXaxis()->SetTitle("mm");
			h1_tpc_x[i][j]->GetYaxis()->SetTitle("count");
			h1_tpc_x[i][j]->SetFillColor(2+i);
			h1_tpc_x[i][j]->SetLineColor(kBlack);
		}
		for(int j=0; j<4; ++j) {
			h1_tpc_y[i][j] = new TH1I(Form("TPC%d_Y%d", i, j), Form("TPC%d y (anode %d)", i, j), 180, -30, 30);
			h1_tpc_y[i][j]->GetXaxis()->SetTitle("mm");
			h1_tpc_y[i][j]->GetYaxis()->SetTitle("count");
			h1_tpc_y[i][j]->SetFillColor(2+i);
			h1_tpc_y[i][j]->SetLineColor(kBlack);
		}
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		for(int i=0; i<N_TPC; ++i) {
			const auto& tpc = frs->tpc[i];
			for(const auto& hit_a4 : tpc.hits) {
				for(int a=0; a<4; ++a) {
					const auto& hit = hit_a4[a];
					int dl = (a >> 1);
					h1_tpc_x[i][dl] -> Fill(hit.x);
					h1_tpc_y[i][a]  -> Fill(hit.y);
				}
			}
		}
	}

	TCanvas *cX[N_TPC];
	TCanvas *cY[N_TPC];
	for(int i=0; i<N_TPC; ++i) {
		cX[i] = new TCanvas(Form("TPC%d-X", i), Form("TPC%d-X", i), 1600, 1200);
		cX[i]->Divide(1,2);

		cX[i]->cd(1);
		h1_tpc_x[i][0]->Draw();

		cX[i]->cd(2);
		h1_tpc_x[i][1]->Draw();

		cY[i] = new TCanvas(Form("TPC%d-Y", i), Form("TPC%d-Y", i), 1600, 1200);
		cY[i]->Divide(2,2);
		
		cY[i]->cd(1);
		h1_tpc_y[i][0]->Draw();
		
		cY[i]->cd(2);
		h1_tpc_y[i][1]->Draw();
		
		cY[i]->cd(3);
		h1_tpc_y[i][2]->Draw();
		
		cY[i]->cd(4);
		h1_tpc_y[i][3]->Draw();
	}
}
