#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

void test_map(std::string fileName = "") {
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSMap>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);

	TH1I* h1 = new TH1I("h1", "h1", 10, 0, 10); 
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		h1->Fill(frs->tpc[0].tdc[0].size());
	}

	h1->Draw();
}
