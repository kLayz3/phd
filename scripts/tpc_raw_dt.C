#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

void tpc_raw_dt(std::string fileName = "", int i=0) {
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH2I* h2_tpc_d = new TH2I(Form("h2_tpc_d%d", i), Form("TPC%d: (x1-x0) vs. (x1+x0)/2", i), 300, -100, 100, 100, -4, 4);
	
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		auto& tpc = frs->tpc.at(i);
		if(tpc.hits.size() != 1) continue;
		
		auto& m = tpc.hits[0];
		if(std::isnan(m[0].x) or std::isnan(m[2].x))
			continue;

		h2_tpc_d->Fill ( 
			(m[2].x + m[0].x) / 2,
			 m[2].x - m[0].x
		);
	}
	h2_tpc_d->GetXaxis()->SetTitle("#frac{x_{1}+x_{0}}{2} [mm]");
	h2_tpc_d->GetYaxis()->SetTitle("x_{1} - x_{0} [mm]");
	TCanvas* c = new TCanvas("c", "c", 2000, 1000);
	h2_tpc_d->Draw("COLZ");

}
