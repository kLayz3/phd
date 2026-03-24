/* Procedure is the following:
 * [1] Gate on very narrow strip around delta = 0 which represents
 * very central hits.
 * [2] On this cut, draw maximum strip value `fCP`, initially binned (640,0,640).
 * This should be even across the detector for one specific ion charge. */

#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/FitSpline.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

constexpr int D_BINS = 1000;
constexpr int N_STRIPS = 64;

enum class DoFit {no, yes, verify};

void foot_gain_match (
	std::string fileName = "", 
	int ifoot = 0,
	int bins_per_asic = 64,
	DoFit do_fit = DoFit::no,
	std::array<double,3> foot_binning = {1000, 4, 4000}, 
	std::array<double,2> delta_cut = {-0.05, 0.05},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX} 
) {
	if(N_STRIPS % bins_per_asic != 0)
		throw std::runtime_error(Form("Passed: %d , not evenly divisible by %d.\n",
			bins_per_asic, N_STRIPS));

	FOOTParam *foot_param; 
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		foot_param = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!foot_param)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
	}
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h1_delta     = new TH1P("((h1_d0))Delta [from -0.5, 0.5]", kGreen-2, 150, -0.499, 0.499); 
	auto* h1_delta_cut = new TH1P("((h1_d1))Delta [from -0.5, 0.5]", kGreen-3, 150, -0.499, 0.499); 
	auto* h1_foot_e = new TH1P("FOOT E [ADC units]@Central strip value", ORGB{0xB2FD30}, (int)(1.5*foot_binning[0]), foot_binning[1], foot_binning[2]); 
	auto* hit_energy = new TH2P(Form("((h2_foot%d))Max Signal [ADC]:Strip number [0..640]@FOOT%d Raw, per ASIC", ifoot, ifoot), 
		bins_per_asic*10, 0,640,
		foot_binning[0], foot_binning[1], foot_binning[2]);

	double C_ADC = FOOTGainParam::CARBON_ADC;
	std::array<double, 3> v_binning = {foot_binning[0], C_ADC/1500 * foot_binning[1], C_ADC/1500 * foot_binning[2]};

	auto* h1_foot_e_corr = new TH1P("((h1_corr))Gain matched FOOT E [Corr ADC units]@Central strip value", ORGB{0x52FD30}, v_binning[0], v_binning[1], v_binning[2]); 
	auto* hit_energy_corr = new TH2P(Form("((h2_foot%d_corr))Max Signal [ADC]:Strip number [0..640]@FOOT%d Corrected, per ASIC", ifoot, ifoot), 
		bins_per_asic*10, 0,640,
		v_binning[0], v_binning[1], v_binning[2]);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];
		if(sci21.hits.size() != 1) continue;
		if(sci22.hits.size() != 1) continue;
		if(sci31.hits.size() != 1) continue;

		h1_sci21->Fill(sci21.E);
		h1_sci22->Fill(sci22.E);
		h1_sci31->Fill(sci31.E);

		if(!mnd::IsInside(sci21.E, sci21_cut)) continue;
		if(!mnd::IsInside(sci22.E, sci22_cut)) continue;
		if(!mnd::IsInside(sci31.E, sci31_cut)) continue;

		h1_sci21_cut->Fill(sci21.E);
		h1_sci22_cut->Fill(sci22.E);
		h1_sci31_cut->Fill(sci31.E);

		for(const auto& cl : foot->fCl) {
			double delta = cl.Delta();
			h1_delta->Fill(delta);
			if(!mnd::IsInside(delta, delta_cut)) continue;
			
			h1_delta_cut->Fill(delta);
			int i = static_cast<int>( cl.fCX );
			double e = cl.fCP;
			hit_energy->Fill(i, e);
			h1_foot_e->Fill(e);
			
			if(do_fit == DoFit::verify) {
				double cog = cl.fCX;
				e /= foot_param->gain.CorrectionFactor(cog);
				hit_energy_corr->Fill(cog, e);
				h1_foot_e_corr->Fill(e);
			}
		}
	}

	std::vector<TGraph*>      profile_fit, gauss_fit;
	std::vector<TGraphErrors*> profile_raw, gauss_raw;

#define POLY_DEG 4
	/* Try to fit a spline(s) for middle 2 ASICs. */
	if(do_fit == DoFit::yes) { 
		//for(int a: {4,5}) {
		for(int a: {2,3,4,5,6,7}) {
		//for(int a=0; a<8; ++a) {
			double x_lo  = (a) * 64 + 0.01;
			double x_hi = (a+1) * 64 - 0.01;
			auto [rg, graw, gfit] = FitSplineAndGraph<POLY_DEG, fit_info::GAUSS_MAX> (
					*hit_energy, x_lo, x_hi, 40, 1.1 /*, Verbosity::CHATTY */
					);
			auto [rp, praw, pfit] = FitSplineAndGraph<POLY_DEG, fit_info::PROFILE_MAX> (
					*hit_energy, x_lo, x_hi, 40 /*, 1.1, Verbosity::CHATTY */
					);
			profile_fit.push_back(pfit);
			profile_raw.push_back(praw);
			gauss_fit.push_back(gfit);
			gauss_raw.push_back(graw);

			std::cout << "FOOT" <<  foot_param->de10_index_ << ": parameters (ASIC: " << a << "): " << rg << std::endl;
		}
	}

	/* Do the small fit in the 1D plot. */
	const double sratio = 0.5;
	auto [fitr, err_] = GaussFitMax(*h1_foot_e, 0.6);
	TH1D* h = *h1_foot_e; 
	printf("1D projection yields: max: %.2f, gauss fit max (around this max+-%.1f sigma): %.2f +- %.2f\n",
		h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ), sratio,
		fitr[1], err_[1]);
	
	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, foot_binning[1], 
				                i * 64, foot_binning[2]);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}
	
	TCanvas *c = new TCanvas(Form("cRAW%d", ifoot), Form("FOOT%d", ifoot), 2000, 1400);
	hit_energy->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1((*hit_energy)->GetYaxis()->GetXmin());
		l->SetY2((*hit_energy)->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	for(auto* pfit : profile_fit) pfit->Draw("L SAME");
	for(auto* praw : profile_raw) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw) graw->Draw("P SAME");

	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* cd = new TCanvas("cd", "Delta", 1400, 800);
	cd->Divide(2,2);
	cd->cd(1); h1_delta->Draw();
	cd->cd(2); h1_delta_cut->Draw();
	cd->cd(3); h1_foot_e->Draw();

	if(do_fit == DoFit::verify) {
		cd->cd(4); h1_foot_e_corr->Draw();
		TCanvas *cc = new TCanvas(Form("cCORR%d", ifoot), Form("[Corr] FOOT%d", ifoot), 2000, 1400);
		TH2D* h2 = &hit_energy_corr->h;
		h2->Draw("COLZ");
		
		for(auto* l0 : vlines) {
			TLine* l = dynamic_cast<TLine*>(l0->Clone());
			l->SetY1(h2->GetYaxis()->GetXmin());
			l->SetY2(h2->GetYaxis()->GetXmax());
			l->Draw("SAME");
		}
	}
}
