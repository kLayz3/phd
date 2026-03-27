#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../../includes/util/PrettyHisto.hxx"
using namespace ROOT::Experimental;
using namespace ROOT;

void pretty_plot_ped (
	std::string fileName = "", 
	int N=20,
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX}
) {
	TFile* f = new TFile(fileName.c_str(), "READ");
	if(!f || f->IsZombie()) { printf("Error opening file.\n"); exit(2); }

#define GET_OBJ(TYPE, x, EXT) \
	auto* x = f->Get<TYPE>(Form("FOOT%d_%s", N, #EXT)); \
	if(!x) { printf("Err getting \'%s\', extension: \'%s\', object: \'%s\', type: \'%s\'\n", \
		#x, #EXT, Form("FOOT%d_%s", N, #EXT), #TYPE); exit(3); } 

	GET_OBJ(TH2I, raw, h2_raw);
	GET_OBJ(TH2D, corr, h2_corr);

	raw->GetXaxis()->SetTitle("Strip number");
	raw->GetYaxis()->SetTitle("ADC value (12-bit)");

	corr->GetXaxis()->SetTitle("Strip number");
	corr->GetYaxis()->SetTitle("Corrected ADC value");
	
	printf("N_STRIPS: %d\nNASIC: %d\n", 64, 10);
	
	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTMap>(Form("FOOT%d", N));
	auto frs = model->MakeField<RNFRSMap>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* corr_cut = new TH2P("((h2_cut))FOOT E [ADC units]:Strip mumber@With sci cut", 640,0,640, 3000, -500, 2500);
	auto* h2_ped_offset = new TH2P("((h2_offset)) Common ASIC offset [ADC units]:ASIC number [0...9]@With cut", 
		10,0,10, 5000, -2500, 2500);
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		if(!foot->HasData()) continue;

		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];
		
		if(sci21.tdc.size() != 1) continue;
		if(sci22.tdc.size() != 1) continue;
		if(sci31.tdc.size() != 1) continue;
		double e21 = sqrt(sci21.qdc[0] * sci21.qdc[1]);	
		double e22 = sqrt(sci22.qdc[0] * sci22.qdc[1]);	
		double e31 = sqrt(sci31.qdc[0] * sci31.qdc[1]);	
		h1_sci21->Fill(e21);
		h1_sci22->Fill(e22);
		h1_sci31->Fill(e31);

		if(!mnd::IsInside(e21, sci21_cut)) continue;
		if(!mnd::IsInside(e22, sci22_cut)) continue;
		if(!mnd::IsInside(e31, sci31_cut)) continue;

		h1_sci21_cut->Fill(e21);
		h1_sci22_cut->Fill(e22);
		h1_sci31_cut->Fill(e31);
		
		for(int i=0; i<640; ++i)
			corr_cut->Fill(i, foot->FOOTE[i]);
		for(int a=0; a<10; ++a) {
			h2_ped_offset->Fill(a, foot->common_offset[a]);
		}
	}

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, -500, 
				                 i * 64, 4096);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}

	TCanvas* cped;
	cped = new TCanvas("cped", "Pedestal", 1000, 2000);
	cped->Divide(3,1);

	cped->cd(1);
	gPad->SetLogz();
	raw->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(raw->GetYaxis()->GetXmin());
		l->SetY2(raw->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}

	cped->cd(2);
	gPad->SetLogz();
	corr->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(corr->GetYaxis()->GetXmin());
		l->SetY2(corr->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	cped->cd(3); gPad->SetLogz();
	corr_cut->Draw("COLZ");

	constexpr int NS = 640;
	typedef std::array<double, NS> Arr; 

	GET_OBJ(Arr, s0_, ped_sigma);
	TGraph* s0 = new TGraph(NS, s0_->data());
	s0->SetMarkerStyle(20);
	s0->SetMarkerSize(1.1);
	s0->SetMarkerColor(kRed - 1);
	TCanvas* csig = new TCanvas("csig", "Sigma", 1618, 1000);
	csig->Divide(2,1);
	csig->cd(1);
	s0->SetTitle(Form("Pedestal Width: FOOT%d;Strip number;ADC value", N));
	s0->Draw("AP");
	csig->cd(2);
	h2_ped_offset->Draw("COLZ");
	
	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

}
