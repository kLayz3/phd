#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../includes/PrettyHisto.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

const char* inv_arg = "Second arg must be either string: H|He|Li|Be|B|C or pair of integers for the range.";

void foot_eta_corr(std::string fileName = "", 
	int ifoot = 0, 
	std::array<double,2> foot_cut = {NAN,NAN}, 
	std::array<double,2> sci21_cut = {NAN,NAN},
	std::array<double,2> sci22_cut = {NAN,NAN} 
) {
	if(isnan(sci21_cut[0])) { sci21_cut[0] = -1; }
	if(isnan(sci21_cut[1])) { sci21_cut[1] = DBL_MAX; }
	if(isnan(sci22_cut[0])) { sci22_cut[0] = -1; }
	if(isnan(sci22_cut[1])) { sci22_cut[1] = DBL_MAX; }

	if(isnan(foot_cut[0]) and isnan(foot_cut[1]) ) {
		std::cout << "Enter the foot_cut arg: H|He|Li|Be|B|C : ";
		
		std::string s{};
		std::cin >> s;

		if(s == "He") {
			foot_cut[0] = 40; foot_cut[1] = 100;
		} else if(s == "Li") {
			throw std::invalid_argument("`Li` not supported yet. It's unclear.");
		} else if(s == "Be") {
			foot_cut[0] = 400; foot_cut[1] = 650;
		} else if(s == "B" ) {
			foot_cut[0] = 650; foot_cut[1] = 1000;
		} else if(s == "C" ) {
			foot_cut[0] = 950; foot_cut[1] = 1450;
		} else {
			throw std::invalid_argument(inv_arg);
		}
		std::cout << "\nOk, taking bounds: " << foot_cut[0] << ", " << foot_cut[1] << endl;
	}
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH2I* h2_eta = new TH2I("h2_eta", Form("Fit Energy (integral) vs. delta FOOT%d", ifoot), 80, -0.55, 0.55, 80, foot_cut[0], foot_cut[1]);
	TH1I* h1_delta = new TH1I("h1_delta", Form("Delta FOOT%d", ifoot), 80, -0.55, 0.55);
	h1_delta->GetYaxis()->SetTitle("count");
	h1_delta->GetXaxis()->SetTitle("Delta (from -0.5 to 0.5)");
	h1_delta->SetFillStyle(1001); h2_eta->SetFillColor(kGreen + 1);

	TH2I* h2_m_vs_delta = new TH2I("h2_m_vs_delta", Form("Multp. vs delta FOOT%d", ifoot), 80, -0.55, 0.55, 10, 0.5, 10.5);
	h2_m_vs_delta->GetYaxis()->SetTitle("Cluster size (integer)");
	h2_m_vs_delta->GetXaxis()->SetTitle("Delta (from -0.5 to 0.5)");

	TH2I* tmp = new TH2I("h2_tmp", "Sigma fit vs. delta", 80,-0.55, 0.55, 100, 0.2, 1.0);
	tmp->GetYaxis()->SetTitle("Sigma of the fitted Gauss");
	tmp->GetXaxis()->SetTitle("Delta (from -0.5 to 0.5)");

	TH2I* sum_minus_fit_v_delta = new TH2I("sum_minus_fit_v_delta", "Cluster sum - Fit integral vs. delta", 80,-0.55, 0.55,
		200, -500, 500);
	sum_minus_fit_v_delta->GetYaxis()->SetTitle("Cluster energy sum - Fitted Gauss Integral (A0*sigma*sqrt(2pi))");
	sum_minus_fit_v_delta->GetXaxis()->SetTitle("Delta (from -0.5 to 0.5)");

	TH2I* sum_energy_vs_delta = new TH2I("sum_energy_vs_delta", "Cluster sum vs. delta", 80,-0.55, 0.55,
		80, foot_cut[0], foot_cut[1]);
	sum_energy_vs_delta->GetYaxis()->SetTitle("Cluster energy sum");
	sum_energy_vs_delta->GetXaxis()->SetTitle("Delta (from -0.5 to 0.5)");

	h2_eta->GetXaxis()->SetTitle("Delta");
	h2_eta->GetYaxis()->SetTitle("E from fit");

	TH1I* h1_delta_diff = new TH1I("h1_delta_diff", Form("Delta(fit) - Delta FOOT%d", ifoot), 150, -0.3, 0.3);
	h1_delta_diff->GetYaxis()->SetTitle("count");
	h1_delta_diff->GetXaxis()->SetTitle("Delta diff (fit-raw)");
	h1_delta_diff->SetFillStyle(1001); h2_eta->SetFillColor(kRed + 1);

	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h2_sci   = new TH2P("SCI21 QDC mean [QDC units]:SCI22 QDC mean [QDC units]", 500, 300, 4000, 500, 300, 4000);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		if(sci21.hits.size() != 1) continue;
		if(sci22.hits.size() != 1) continue;
		
		(*h1_sci21)->Fill(sci21.E);
		(*h1_sci22)->Fill(sci22.E);
		h2_sci->Fill(sci21.E, sci22.E);

		if(sci21.E < sci21_cut[0] or sci21.E > sci21_cut[1]) continue;
		if(sci22.E < sci22_cut[0] or sci22.E > sci22_cut[1]) continue;

		for(const auto& cl : foot->fCl) {
			if(!cl.fit.IsOk()) continue;

			double e = cl.fit.E();
			double mean_position = cl.fCX;
			double d0 = (double)mnd::rround<int>(mean_position) - mean_position;

			if(foot_cut[0] <= e and e < foot_cut[1]) {
				h2_eta->Fill( d0, e );
				h1_delta->Fill( d0 );
				h2_m_vs_delta->Fill(d0 , cl.fCM);
				tmp->Fill(d0, cl.fit.sigma);
				sum_minus_fit_v_delta -> Fill(d0, cl.fCE - cl.fit.E() );
				sum_energy_vs_delta -> Fill(d0, cl.fCE);

				h1_delta_diff->Fill( cl.fit.delta - d0 );
			}
		}
	}

	TCanvas *c = new TCanvas(Form("cRAW%d", ifoot), Form("eta%d", ifoot), 2000, 1400);
	c->Divide(3,2);
	c->cd(1); gPad->SetGrid();
	h2_eta->Draw("COLZ");
	c->cd(2); gPad->SetGrid();
	h1_delta->Draw();
	c->cd(3); gPad->SetGrid();
	h2_m_vs_delta->Draw("COLZ");

	c->cd(4); gPad->SetGrid();
	tmp->Draw("COLZ");
	c->cd(5); gPad->SetGrid();
	sum_minus_fit_v_delta->Draw("COLZ");
	c->cd(6); gPad->SetGrid();
	sum_energy_vs_delta->Draw("COLZ");
	
	//TCanvas* cdiff = new TCanvas(Form("cddiff%d", ifoot), Form("Diff delta(fit) - delta FOOT%d", ifoot), 1200, 600);
	//h1_delta_diff->Draw(); gPad->SetGrid();

	TCanvas* cs = new TCanvas("cs", "cs", 1600, 1200);
	cs->Divide(1,2);
	cs->cd(1); (*h1_sci21)->Draw(); gPad->SetGrid();
	cs->cd(2); (*h1_sci22)->Draw(); gPad->SetGrid();

	TCanvas* cs2 = new TCanvas("cs2", "cs2", 1600, 1200);
	h2_sci->Draw("COLZ"); gPad->SetGrid();
}
