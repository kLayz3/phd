#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

constexpr int xs[] = {1,3,5,7};
constexpr int ys[] = {0,2,4,6};
constexpr size_t NFOOT = 8;
const char* orientation[] = {
	"y", "x", "y", "x", "y", "x", "y", "x" 
};

void foot_fit_explo(std::string fileName = "", int ifoot = -1) {
	if(ifoot != -1 and ifoot > NFOOT)
		throw std::invalid_argument(Form("Second arg must be <%zu to plot specific FOOT.\n", NFOOT));

	ROOT::EnableImplicitMT();

	TH1I *h1_foot_x[NFOOT];
	TH1I *h1_foot_eC[NFOOT];
	TH1I *h1_foot_eD[NFOOT];
	TH1I *h1_foot_s[NFOOT];
	TH1I* h1_foot_xnofit[NFOOT];
	TH1I* h1_foot_enofit[NFOOT];
	TH1I* h1_foot_efull[NFOOT];
	TH2I* h2_foot_mult[NFOOT];

	for(int i=0; i<NFOOT; ++i) {
		h1_foot_x[i] = new TH1I(Form("FOOT%d_x", i), Form("FOOT%d measuring \'%s\'", i, orientation[i]), 640,0,640);
		h1_foot_x[i]->GetXaxis()->SetTitle(Form("%s in units of strip#", orientation[i]));
		h1_foot_x[i]->GetYaxis()->SetTitle("Count");
		h1_foot_x[i]->SetFillStyle(1001); h1_foot_x[i]->SetFillColor(kGreen + 1);

		h1_foot_eC[i] = new TH1I(Form("FOOT%d_eC", i), Form("FOOT%d energy (continuous extrapolation)", i), 2000,0,2000);
		h1_foot_eC[i]->GetXaxis()->SetTitle("So-called \'energy\' (ADC units sum)");
		h1_foot_eC[i]->GetYaxis()->SetTitle("Count");
		h1_foot_eC[i]->SetFillStyle(1001); h1_foot_eC[i]->SetFillColor(kYellow - 3);

		h1_foot_eD[i] = new TH1I(Form("FOOT%d_eD", i), Form("FOOT%d energy (discrete extrapolation)", i), 2000,0,2000);
		h1_foot_eD[i]->GetXaxis()->SetTitle("So-called \'energy\' (ADC units sum)");
		h1_foot_eD[i]->GetYaxis()->SetTitle("Count");
		h1_foot_eD[i]->SetFillStyle(1001); h1_foot_eD[i]->SetFillColor(kYellow - 2);

		h1_foot_s[i] = new TH1I(Form("FOOT%d_s", i), Form("FOOT%d cluster fit width", i), 200, 0, 3);
		h1_foot_s[i]->GetXaxis()->SetTitle("Width in units of strips");
		h1_foot_s[i]->GetYaxis()->SetTitle("Count");
		h1_foot_s[i]->SetFillStyle(1001); h1_foot_s[i]->SetFillColor(kCyan - 2);

		h1_foot_xnofit[i] = new TH1I(Form("FOOT%d_xnofit", i), Form("FOOT%d measuring \'%s\' (cluster size=1,2)", i, orientation[i]), 640,0,640);
		h1_foot_xnofit[i]->GetXaxis()->SetTitle(Form("%s in units of strip#", orientation[i]));
		h1_foot_xnofit[i]->GetYaxis()->SetTitle("Count");
		h1_foot_xnofit[i]->SetFillStyle(1001); h1_foot_xnofit[i]->SetFillColor(kGreen - 2);

		h1_foot_enofit[i] = new TH1I(Form("FOOT%d_enofit", i), Form("FOOT%d energy (cluster size=1,2)", i),600,0,200);
		h1_foot_enofit[i]->GetXaxis()->SetTitle("Energy (ADC units)");
		h1_foot_enofit[i]->GetYaxis()->SetTitle("Count");
		h1_foot_enofit[i]->SetFillStyle(1001); h1_foot_enofit[i]->SetFillColor(kRed - 2);

		h1_foot_efull[i] = new TH1I(Form("FOOT%d_efull", i), Form("FOOT%d total energy", i),6000,0,2000);
		h1_foot_efull[i]->GetXaxis()->SetTitle("Energy (ADC units)");
		h1_foot_efull[i]->GetYaxis()->SetTitle("Count");
		h1_foot_efull[i]->SetFillStyle(1001); h1_foot_efull[i]->SetFillColor(kRed + 2);

		h2_foot_mult[i] = new TH2I(Form("FOOT%d_mult", i), Form("FOOT%d multiplicity vs. cluster energy (cont.extr.)", i),
				1000,0,2000, 10, 0.5, 10.5);
		h2_foot_mult[i]->GetXaxis()->SetTitle("Cluster energy derived from continuous extrapolation");
		h2_foot_mult[i]->GetYaxis()->SetTitle("Absolute multiplicity");
	}

	/* Quickly also fetch objects from file. 
	std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
	for(int i=0; i<NFOOT; ++i) {
		h2_foot_mult[i] = f->Get<TH2I>(Form("FOOT%d_h2_mult_e", i));
		h2_foot_mult[i]->SetDirectory(nullptr);
		h2_foot_mult[i]->RebinX(32);
	} f.reset(nullptr); */

	auto model = RNTupleModel::Create();
	std::array <
		std::shared_ptr<RNFOOTCal>, NFOOT
	> foot {};
	
	for(int i=0; i<NFOOT; ++i) {
		foot[i] = model->MakeField<RNFOOTCal>(Form("FOOT%d", i));
	}
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		for(int i=0; i<NFOOT; ++i) {
			if(ifoot >= 0 and i != ifoot) continue;
			
			for(const auto& cl : foot[i]->fCl) {
				const FOOTClusterFit& cl_fit = cl.fit;
					
				if(cl_fit.IsOk()) {
					double e = cl_fit.E();
					h1_foot_x[i] ->Fill( cl_fit.X() );
					h1_foot_eC[i]->Fill( e );
					h1_foot_eD[i]->Fill(cl_fit.E_discrete() );
					h1_foot_s[i] ->Fill( cl_fit.sigma );
					h2_foot_mult[i]->Fill(e , cl.fCM);
				}
				else {
					h1_foot_xnofit[i]->Fill( cl.fCX ); 
					h1_foot_enofit[i]->Fill( cl.fCE ); 
					h2_foot_mult[i]->Fill( cl.fCE , cl.fCM);
				}
				h1_foot_efull[i]->Fill( cl.fCE ); 
			}
		}
	}
	
	
	for(int i=0; i<NFOOT; ++i) {
		if(ifoot >= 0 and i != ifoot) continue;
		
		TCanvas *c = new TCanvas(Form("c%d", i), Form("c%d", i), 1900, 1250);
		c->Divide(4,2);
#define HILFE_(h,i, ...) \
		c->cd(i); \
		gPad->SetGrid(); \
		h->Draw(__VA_ARGS__);
		
		HILFE_(h1_foot_eC[i], 1, "HIST")
		HILFE_(h1_foot_eD[i], 2, "HIST")
		HILFE_(h1_foot_s[i], 4, "HIST")
		HILFE_(h1_foot_x[i], 5, "HIST")
		HILFE_(h1_foot_xnofit[i], 6, "HIST")
		HILFE_(h1_foot_enofit[i], 3, "HIST")
		HILFE_(h1_foot_efull[i], 7, "HIST")
		HILFE_(h2_foot_mult[i], 8, "COLZ")
	
		c->cd(8);
		// Horizontal dashed line
		TLine *line = new TLine(0, 2.5, 2000, 2.5);
		line->SetLineStyle(2);
		line->SetLineWidth(6);
		line->SetLineColor(kRed);
		line->Draw("SAME");
		gPad->SetLogz();
	}

	if(ifoot >= 0) {
		auto fname_short = std::string(strrchr(fileName.c_str(), '/') ? strrchr(fileName.c_str(), '/') + 1 : fileName);
		TCanvas* c = (TCanvas*)gPad->GetCanvas();
		if(!c) fprintf(stderr, "gpad static cast to tcanvas failed?"); 
		c->SetTitle(fname_short.substr(0, fname_short.find(".root")).c_str());
		c->Modified();
		c->Update();
	}
}
