void pretty_plot_ped(const char* fileName = "", int N=20) {
	TFile* f = new TFile(fileName, "READ");
	if(!f || f->IsZombie()) { printf("Error opening file.\n"); exit(2); }

#define GET_OBJ(TYPE, x, EXT) \
	auto* x = dynamic_cast<TYPE*>(f->Get(Form("FOOT%d_%s", N, #EXT))); \
	if(!x) { printf("Err getting \'%s\', extension: \'%s\', object: \'%s\', type: \'%s\'\n", \
		#x, #EXT, Form("FOOT%d_%s", N, #EXT), #TYPE); exit(3); } 

	GET_OBJ(TH2I, raw, h2_raw);
	GET_OBJ(TH2D, mid, h2_mid);
	GET_OBJ(TH2D, corr, h2_corr);

	raw->GetXaxis()->SetTitle("Strip number");
	raw->GetYaxis()->SetTitle("ADC value (12-bit)");

	mid->GetXaxis()->SetTitle("Strip number");
	mid->GetYaxis()->SetTitle("ADC value minus global pedestal");

	corr->GetXaxis()->SetTitle("Strip number");
	corr->GetYaxis()->SetTitle("Corrected ADC value");
	
	printf("N_STRIPS: %d\nNASIC: %d\n", 64, 10);

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
	
	cped->Divide(1,3);

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
	mid->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(mid->GetYaxis()->GetXmin());
		l->SetY2(mid->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}

	cped->cd(3);
	gPad->SetLogz();
	corr->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(corr->GetYaxis()->GetXmin());
		l->SetY2(corr->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}


	GET_OBJ(TGraph, s0, sigma0);
	GET_OBJ(TGraph, s1, sigma1);

	TCanvas* csig;
	csig = new TCanvas("csig", "Sigma", 1618, 1000);

	csig->cd();
	s0->SetTitle(Form("Pedestal Width: FOOT%d;Strip number;ADC value", N));
	s0->Draw("AP");
	s1->Draw("P SAME");
	
	double _ymin = std::min(
			TMath::MinElement(s0->GetN(), s0->GetY()),
			TMath::MinElement(s1->GetN(), s1->GetY())
		);
	double _ymax = std::max(
			TMath::MaxElement(s0->GetN(), s0->GetY()),
			TMath::MaxElement(s1->GetN(), s1->GetY())
		);
	double _range = _ymax - _ymin;
	printf("Min: %.2f; max: %.2f; range: %.2f\n", _ymin, _ymax, _range);
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(_ymin - 0.09 * _range);
		l->SetY2(_ymax - 0.09 * _range);
		l->SetLineColor(kBlack);
		l->SetLineWidth(2);
		l->Draw("SAME");
	}

	TLegend *leg = new TLegend(0.77, 0.77, 0.9, 0.9);
	leg->AddEntry(s0, "Without Fine Sub", "p");
	leg->AddEntry(s1, "With Fine Sub", "p");
	leg->SetTextFont(42);
	leg->SetTextSize(0.022);
	leg->SetLineWidth(2);
	leg->Draw();
}
