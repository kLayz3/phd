#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

void foot_spread(std::string fileName = "") {
	ROOT::EnableImplicitMT();
	
	auto model = RNTupleModel::Create();
	std::array <
		std::shared_ptr<RNFOOTCal>, 8
	> foot {};
	for(int i=0; i<8; ++i) {
		foot[i] = model->MakeField<RNFOOTCal>(Form("FOOT%d", i));
	}
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH1I *h2_foot[8];
	for(int i=0; i<8; ++i) {
		h2_foot[i] = new TH1I(Form("FOOT%d_x", i), Form("FOOT%d_x", i), 150,0,640);
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		for(int i=0; i<8; ++i) {
			std::vector<double> cl_pos = foot[i]->X();
			for(const auto p : cl_pos) {
					h2_foot[i]->Fill(p);
			}
		}
	}
	
	/* Try to fit a Gauss around them. */
	TCanvas* c[8];
	for(int i=0; i<8; ++i) {
		TH1I* foot = h2_foot[i];
		int maxBin = foot->GetMaximumBin();
		double v = foot->GetBinContent(maxBin);
	
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 1600, 1200);
	
		h2_foot[i]->Draw("COLZ");
		//foot->Fit("gaus", "Q", "", 
	}



}
