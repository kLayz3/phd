#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"

using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;

void pretty_cal(const char* fileName = "", int N = 4, uint32_t Nevents = 16, uint64_t nFirstEntry = 0) {
	TCanvas *cbad = new TCanvas("cbad", "Bad events", 2600, 2000);

	const int N_DIV = 2;
	const int N_DIV2 = N_DIV * N_DIV; 
	const int N_GRAPHS_PER_PAD = 1; /* 1,2,3,4. */
	Nevents = (N_DIV2 * N_GRAPHS_PER_PAD < 64) ? (N_DIV2 * N_GRAPHS_PER_PAD) : 64;

	cbad->Divide(N_DIV, N_DIV);
	
	int curr_ev = 0;

	int style[]{22, 23, 21, 47};
	int col[]{2, 3, 4, 6};

	double arr[640]{};
	std::iota(arr, arr+640,0);
	
	auto model = RNTupleModel::Create();
	auto footc = model->MakeField<RNFOOTCal>(Form("FOOT%d", N)); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
	uint64_t nentries = ntuple->GetNEntries();
	for(uint64_t entryId = nFirstEntry; entryId < nentries; ++entryId) {
		ntuple->LoadEntry(entryId);
		if(footc->_fBadE.size() == 0) continue;
		
		if(entryId == Nevents) break;

		TGraph* gr = new TGraph(640);
		memcpy(gr->GetX(), arr, sizeof(arr));
		memcpy(gr->GetY(), footc->_fBadE.data(), sizeof(arr));
		
		gr->SetMarkerStyle(style[curr_ev % N_GRAPHS_PER_PAD]);
		gr->SetMarkerColor(col[curr_ev % N_GRAPHS_PER_PAD]);
		gr->SetMarkerSize(0.55);
		
		if(curr_ev % N_GRAPHS_PER_PAD == 0) {
			cbad->cd(curr_ev / N_GRAPHS_PER_PAD + 1);
			gr->SetTitle(Form("FOOT%d;StripNum;ADC(ca)", N));
			gr->Draw("AP");
			TLatex *latex = new TLatex{};
			latex->SetTextAlign(13);  //align at top
			latex->DrawLatex(520, 40, Form("Ev: %lu", entryId));  
		} else
			gr->Draw("P SAME");
		++curr_ev;
	}

	TCanvas *cweird = new TCanvas("cweird", "Weird events", 2600, 2000);
	cweird->Divide(N_DIV, N_DIV);

	curr_ev = 0;

	srand(time(NULL));
	for(uint64_t entryId = nFirstEntry; entryId < nentries; ++entryId) {
		ntuple->LoadEntry(entryId);
		if(footc->_fHeClSize1.size() == 0) continue;
		if(curr_ev == Nevents) break;

		TGraph* gr = new TGraph(640);
		memcpy(gr->GetX(), arr, sizeof(arr));
		memcpy(gr->GetY(), footc->_fHeClSize1.data(), sizeof(arr));
		
		gr->SetMarkerStyle(style[curr_ev % N_GRAPHS_PER_PAD]);
		gr->SetMarkerColor(col[curr_ev % N_GRAPHS_PER_PAD]);
		gr->SetMarkerSize(0.75);
			
		if(curr_ev % N_GRAPHS_PER_PAD == 0) {
			cweird->cd(curr_ev / N_GRAPHS_PER_PAD + 1);
			gr->SetTitle(Form("FOOT%d;StripNum;ADC(cal)", N));
			gr->Draw("AP");
			TLatex *latex = new TLatex{};
			latex->SetTextAlign(13);  //align at top
			latex->DrawLatex(520, 200, Form("Ev: %lu", entryId));  
		}
		TGraph *gr_cl = new TGraph(footc->fCl.size());
		memcpy(gr_cl->GetX(), footc->X().data(), N*sizeof(double));
		memcpy(gr_cl->GetY(), footc->E().data(), N*sizeof(double));

		gr_cl->Draw("P SAME");
		gr_cl->SetMarkerStyle(29);
		gr_cl->SetMarkerColor(6);
		gr_cl->SetMarkerSize(0.8);

		++curr_ev;
	}
	ntuple.reset();
	footc.reset();

	// leak it intentionally, otherwise objects are dead upon termination.
	TFile* f = TFile::Open(fileName, "READ"); 
	if(!f || f->IsZombie()) { printf("Error opening file.\n"); exit(2); }

	TCanvas *csn = new TCanvas("csn", "csn", 2600,2000);
	csn->Divide(2, 4);

#define GET_OBJ(TYPE, x, N, EXT) \
	x = f->Get<TYPE>(Form("FOOT%d_%s", N, #EXT)); \
	if(!x) { printf("Err getting \'%s\', extension: \'%s\', object: \'%s\', type: \'%s\'\n", \
		#x, #EXT, Form("FOOT%d_%s", N, #EXT), #TYPE); exit(3); }
	TH1I *hsn[8];
	GET_OBJ(TH1I, hsn[0], 0, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[1], 1, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[2], 2, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[3], 3, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[4], 4, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[5], 5, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[6], 6, h1_sn_ratio);
	GET_OBJ(TH1I, hsn[7], 7, h1_sn_ratio);

	for(int i=0; i<8; ++i) {
		csn->cd(i+1);
		auto* h = hsn[i];
		h->Draw();
		h->GetXaxis()->SetTitle("#sigma_{<}*#sigma_{>}/#sigma_0^2");
		h->GetYaxis()->SetTitle("Count");
		double mx = h->GetBinCenter( h->GetMaximumBin() );
		double my = h->GetBinContent( h->GetMaximumBin() );
		TLine *l = new TLine(mx, 0, mx, my);
		l->SetLineColor(kRed);
		l->SetLineStyle(2);
		l->SetLineWidth(3);
		l->Draw("SAME");
		TLatex *latex = new TLatex{};
		latex->SetTextAlign(13);  //align at top
		latex->DrawLatex(mx + 0.4, my, Form("%.3f", mx)); 
	}
}
