#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "../../includes/PrettyHisto.hxx"

using namespace ROOT::Experimental;

//static const char* label[] = {
//	"21", "22", "23", "24", "41", "42", "31"
//};


void tpc_limits_exploration(const char* fileName = "", int i_tpc=0) {
	if(i_tpc >= RNFRSCal::N_VALID_TPC)
		throw std::invalid_argument(Form("TPC index must be < %d\n", RNFRSCal::N_VALID_TPC));

	const auto& label = RNFRSCal::tpc_label;
	ROOT::EnableImplicitMT();

	TH1P* h1_dl_left_diff_lim = new TH1P(Form("DL1 - DL0 left [ADC units]@TPC%s", label[i_tpc]), ORGB{0x0070DD}, 2000, -2000,2000);
	TH1P* h1_dl_right_diff_lim = new TH1P(Form("DL1 - DL0 right [ADC units]@TPC%s", label[i_tpc]), ORGB{0xAB2D2D}, 2000, -2000,2000);
	TH1P* h1_csum[4] ;
	TH1P* h1_sci_ref_lim[4];

	for(int a=0; a<4; ++a) {
		h1_csum[a] = new TH1P(Form("Anode %d CSUM [ADC units]@TPC%s", a, label[i_tpc]), ORGB{0xFC30FC}, 30000, 0, 30000);
		h1_sci_ref_lim[a] = new TH1P(Form("Anode %d - ref. SCI [ADC units]@TPC%s", a, label[i_tpc]), ORGB{0x03E7FF}, 40000, 0, 40000);
	}

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSMap>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		const auto& tpc = frs->tpc[i_tpc];
		if(tpc.tdc_ref.size() != 1) continue;
		
		i32 sci_ref = tpc.tdc_ref[0];

		for(int d: {0,1}) {
			const auto& hits = tpc.tdc[d];
			if(hits.size() != 1) continue;
			const RNTPCMap::Measurement& hit = hits[0];

			if( std::isfinite(hit.tdc_l) and
				std::isfinite(hit.tdc_r) and 
				std::isfinite(hit.tdc_a[0]) and 
				std::isfinite(hit.tdc_a[1])
			) {
				for(int a: {0,1}) {
					h1_sci_ref_lim[2*d + a] -> Fill( hit.tdc_a[a] - sci_ref );	
					h1_csum[2*d + a] -> Fill( hit.tdc_l + hit.tdc_r - 2 * hit.tdc_a[a] );
				}
			}
		}
		if(tpc.tdc[0].size() == 1 and tpc.tdc[1].size() == 1) {
			h1_dl_left_diff_lim  -> Fill( tpc.tdc[0][0].tdc_l -  tpc.tdc[1][0].tdc_l );
			h1_dl_right_diff_lim -> Fill( tpc.tdc[0][0].tdc_r -  tpc.tdc[1][0].tdc_r );
		}
	}

	TCanvas* c= new TCanvas(Form("c%s", label[i_tpc]), Form("c%s", label[i_tpc]), 1800, 1200);
	c->Divide(4,2);
	for(int a=0; a<4; ++a) {
		c->cd(1+a);
		h1_csum[a]->Draw(); gPad->SetLogy();
		c->cd(5+a);
		h1_sci_ref_lim[a]->Draw(); gPad->SetLogy();
	}

	TCanvas* clr = new TCanvas("clr", "clr", 1800,1200);
	clr->Divide(2,1);
	clr->cd(1); h1_dl_left_diff_lim ->Draw(); gPad->SetLogy();
	clr->cd(2); h1_dl_right_diff_lim->Draw(); gPad->SetLogy();
}
