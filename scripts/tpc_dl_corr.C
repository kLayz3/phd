#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

void tpc_dl_corr(std::string fileName = "", int i=0) {
	ROOT::EnableImplicitMT();
	
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH2I *h1_tpc_xx[2], *h1_tpc_yy[4];
	for(int j=0; j<2; ++j) {
		h1_tpc_xx[j] = new TH2I(Form("TPC%d_X%d_corr", i,j), Form("TPC%d X-correlation (dline %d)", i, j), 
				600, -6, 6, 800, -40, 40);
		h1_tpc_xx[j]->GetXaxis()->SetTitle("X0-X1 [mm]"); 
		h1_tpc_xx[j]->GetYaxis()->SetTitle(Form("X%d [mm]",j)); 
	}
	for(int j=0; j<4; ++j) {
		h1_tpc_yy[j] = new TH2I(Form("TPC%d_Y%d_corr", i, j), Form("TPC%d Y-correlation (anode %d vs. avg.)", i, j), 
				1000, -10, 10, 1000, -40, 40);
		h1_tpc_yy[j]->GetXaxis()->SetTitle(Form("Y%d - Yavg [mm]", j)); 
		h1_tpc_yy[j]->GetYaxis()->SetTitle(Form("Y%d [mm]", j)); 
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		const auto& tpc = frs->tpc[i];
		if(tpc.hits.size() != 1) continue;
		auto& m = tpc.hits[0];
		if(!std::isnan(m[0].x) and !std::isnan(m[2].x)) {
			h1_tpc_xx[0]->Fill(m[0].x - m[2].x, m[0].x);
			h1_tpc_xx[1]->Fill(m[0].x - m[2].x, m[2].x);
		}
		if( !std::isnan(m[0].y) and 
			!std::isnan(m[2].y) and
			!std::isnan(m[2].y) and
			!std::isnan(m[2].y) ) {
			
			double y0 = ( m[0].y + m[1].y + m[2].y + m[3].y ) / 4;
			for(int j=0; j<4; ++j) {
				h1_tpc_yy[j]->Fill(m[j].y - y0, m[j].y);
			}
		}
	}

	TCanvas* cx = new TCanvas("cx", "cx", 2000, 1400);
	cx->Divide(2,1);
	for(int j=0; j<2; ++j) {
		cx->cd(j+1); h1_tpc_xx[j]->Draw("COLZ");
	}

	TCanvas* cy = new TCanvas("cy", "cy", 2000, 1400);
	cy->Divide(2,2);
	for(int j=0; j<4; ++j) {
		cy->cd(j+1); h1_tpc_yy[j]->Draw("COLZ");
	}
}
