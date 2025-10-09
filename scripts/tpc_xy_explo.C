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

	constexpr auto N_TPC = RNFRSCal::N_VALID_TPC;
	TH2I* h1_tpc_x[N_TPC];
	TH2I* h1_tpc_y[N_TPC];

	for(int i=0; i<N_TPC; ++i) {
		h1_tpc_x[i] = new TH2I(Form("TPC%d_X", i), Form("TPC%d x", i), 2,0,2, 200, -100, 100);
		h1_tpc_y[i] = new TH2I(Form("TPC%d_Y", i), Form("TPC%d y", i), 4,0,4, 200, -100, 100);
	}
	uint64_t nvalid = 0;
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		printf("Caught %lu hit\n", ++nvalid);
		
		for(int i=0; i<N_TPC; ++i) {
			const auto& tpc = frs->tpc[i];
			for(const auto& p : tpc.val) {
				printf("Caught TPC%d - %lu hit\n", i, ++nvalid);
				assert(p.delay_line_i < 2 && Form("%d delay_line_i >= 2 ?? Is: %d", i, p.delay_line_i));
				assert(p.anode_i < 4      && Form("%d anode_i >= 4 ?? Is: %d"     , i, p.anode_i));

				h1_tpc_x[i]->Fill(p.delay_line_i, p.x);
				h1_tpc_x[i]->Fill(p.anode_i, p.y);
			}
		}
	}

	TCanvas *c[N_TPC];
	for(int i=0; i<N_TPC; ++i) {
		c[i] = new TCanvas(Form("TPC-%d", i), Form("TPC-%d", i), 1600, 1200);
		c[i]->Divide(1,2);

		c[i]->cd(1);
		h1_tpc_x[i]->Draw("COLZ");

		c[i]->cd(2);
		h1_tpc_y[i]->Draw("COLZ");
	}
}
