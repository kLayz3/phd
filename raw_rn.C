#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TError.h"

using namespace ROOT;
using namespace ROOT::Experimental;

void raw_rn(std::string fileName = "") {
	//gErrorIgnoreLevel = kSysError;
	using TRaw = std::vector<unsigned char>;

	auto model = RNTupleModel::Create();
	auto raw_p = model->MakeField<TRaw>("FRSSortEvent");
	auto ntuple = RNTupleReader::Open(std::move(model), "SortxTree", "./rn/rn_main_0120_sort.root");

	uint64_t nentries = ntuple->GetNEntries();
	printf("Found: %lu entries\n", nentries);

	TH2I *h2 = new TH2I("h2", "h2", 640,0,640, 3000,0,3000);

	for(uint64_t entryId = 0; entryId < nentries; ++entryId) {
		ntuple->LoadEntry(entryId);
		TBufferFile buf(TBuffer::kRead, raw_p->size(),
			raw_p->data(), kFALSE);

		TFRSSortEvent* frs = static_cast<TFRSSortEvent*>( buf.ReadObject(TFRSSortEvent::Class()) );
		
		if(frs->FOOT20 == 640) {
			for(int i=0; i<640; ++i) {
				h2->Fill( frs->FOOT20I[i], frs->FOOT20E[i] );
			}
		}
		delete frs;
	}

	h2->Draw("COLZ");
}
