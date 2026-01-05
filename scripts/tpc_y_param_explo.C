#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"

using namespace ROOT::Experimental;

void tpc_y_param_explo(const char* fileName = "", int i=0) {
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSMap>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);

	TH1I* h1_tpc_xdiff[2];

	for(int a=0; a<2; ++a) {
		h1_tpc_xdiff[a] = new TH1I(Form("TPC%d_XDiff%d", i,a), Form("TPC%d dl(%d) l-r (mult == 1 in both)", i, a), 40000, -20000, 20000);
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		auto& tpc = frs->tpc[i];
		if(tpc.tdc.size() != 1) continue;

		auto& hit = tpc.tdc[0]; // TDC::Measurement;
		for(int j=0; j<2; ++j) {
			if( !std::isnan(hit.tdc_l[j]) and
				!std::isnan(hit.tdc_r[j]) ) 
			{
				h1_tpc_xdiff[j]->Fill (
					hit.tdc_l[j] - hit.tdc_r[j]
				);
			}
		}
	}

	TCanvas* c= new TCanvas(Form("c%d", i), Form("c%d", i), 2000, 1000);
	c->Divide(2,1);
	for(int a=0; a<2; ++a) {
		c->cd(a+1);
		h1_tpc_xdiff[a]->Draw();
	}
}
