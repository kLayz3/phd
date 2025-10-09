#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"

using namespace ROOT::Experimental;

void tpc_y_param_explo(const char* fileName = "") {
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSMap>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);

	constexpr int N_TPC = 4;
	TH1I* h1_tpc_ydiff[N_TPC][4];

	for(int i=0; i<N_TPC; ++i)
		for(int a=0; a<4; ++a)
			h1_tpc_ydiff[i][a] = new TH1I(Form("TPC%d_YDiff%d", i,a), Form("TPC%d anode(%d) - ref (mult == 1 in both)", i, a), 5000, 0, 40000);

	TH1I* clarity_tdc_diff = new TH1I(Form("TDC (n) - (n-1) hit"), Form("TDC (n) - (n-1) hit"), 2000, -100000, 100000);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		for(int i=0; i<N_TPC; ++i) {
			auto& tpc = frs->tpc[i];
			if(tpc.tdc_ref.size() != 1) continue;

			for(int a=0; a<4; ++a)
				if(tpc.tdc_a[a].size() == 1)
					h1_tpc_ydiff[i][a]->Fill(
						tpc.tdc_a[a][0] -
						tpc.tdc_ref[0]
					);
		}

		
		for(int i=0; i<N_TPC; ++i) {
			auto& tpc = frs->tpc[i];
#define TRY_CHN(X) \
			for(const auto& v : tpc.X) { \
				for(int j=1; j < (int)v.size(); ++j) \
					clarity_tdc_diff->Fill(v[j] - v[j-1]); \
			} \
			
			TRY_CHN(tdc_l)
			TRY_CHN(tdc_r)
			TRY_CHN(tdc_a)
		}
	}

	TCanvas* c[N_TPC];
	for(int i=0; i<N_TPC; ++i) {
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 2000, 1000);
		c[i]->Divide(2,2);
		for(int a=0; a<4; ++a) {
			c[i]->cd(a+1);
			h1_tpc_ydiff[i][a]->Draw();
		}
	}

	TCanvas* idk = new TCanvas("idk", "idk", 1600, 1000);
	clarity_tdc_diff->Draw();
}
