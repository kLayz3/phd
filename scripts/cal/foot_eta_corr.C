#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/GaussFitMax.hxx"

#if 0
#	define USING_POLY
#else
#	define USING_FFT
#endif

#ifdef USING_FFT
#	define __USING_LUSTRE_HPC__
#	ifdef __USING_LUSTRE_HPC__
		R__ADD_LIBRARY_PATH(/u/mbajzek/.local/lib)
		R__LOAD_LIBRARY(libfftw3.so)
#	endif
#	include "../../includes/util/FFT.h"
#else
#	include "../../includes/util/FitSpline.hxx"
#	include "../../includes/util/PolyFitter.h"
#endif

using namespace ROOT;
using namespace ROOT::Experimental;

enum class DoGainMatch {no,yes};

/* Because cling issues the WEIRDEST compiler error,..
 * I cannot just simply return `FOOTDeltaParam` instance. It tries to compile the class
 * from the inputs but loses it on template spec in boost preproc library LOL. PEGI 18. */
class FOOTDeltaParam;
FOOTDeltaParam* GetDeltaParams(TH2D*);

constexpr double D_LO = 0.4999;
constexpr int D_BINS = 60;

void foot_eta_corr (
	std::string fileName = "", 
	int ifoot = 0,
	std::array<double,2> foot_cut = {4,4000}, 
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX} 
) {
#ifndef __USING_LUSTRE_HPC__
	gSystem->Load ("libfftw3.so");
#endif

	FOOTParam* p;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		p = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!p)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
	}
	const double CA = FOOTGainParam::CARBON_ADC;

	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH1P* h1_delta = new TH1P(Form("((h1_foot%d))Delta@FOOT%d", ifoot, ifoot), kGreen+1, D_BINS, -D_LO, D_LO);
	TH2P* h2_m_vs_delta = new TH2P(Form("((h2_footm%d))Partial cluster size:delta@FOOT%d", ifoot, ifoot), D_BINS, -D_LO, D_LO, 10, 0.5, 10.5);
	TH2P* sum_energy_vs_delta = new TH2P(Form("((h2_footraw%d))Cluster sum [ADC]:Delta@FOOT%d raw", ifoot, ifoot), 
		D_BINS,-D_LO, D_LO,
		D_BINS, foot_cut[0], foot_cut[1]);
	TH2P* corr_energy_vs_delta = new TH2P(Form("((h2_footcorr%d))E1/%.1f - 1:Delta@FOOT%d", ifoot, CA, ifoot), 
		D_BINS,-D_LO, D_LO,
		D_BINS, foot_cut[0]/CA - 1, foot_cut[1]/CA - 1);
	
	TH1P* h1_m0 = new TH1P(Form("((h_0)) FOOT%d dE [ADC units]@All cluster sizes", ifoot), kYellow - 7, D_BINS, foot_cut[0], foot_cut[1]);
	TH1P* h1_m1 = new TH1P(Form("((h_1)) FOOT%d dE [ADC units]@Only cluster size 1", ifoot), kYellow + 1, D_BINS, foot_cut[0], foot_cut[1]);
	TH1P* h1_m2 = new TH1P(Form("((h_2)) FOOT%d dE [ADC units]@Only cluster sizes >2", ifoot), kYellow - 3, D_BINS, foot_cut[0], foot_cut[1]);

	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_sci    = new TH2P("SCI22 QDC mean [QDC units]:SCI21 QDC mean [QDC units]", 500, 300, 4000, 500, 300, 4000);
	auto* h2_sci3v2 = new TH2P("SCI31 QDC mean [QDC units]:SCI22 QDC mean [QDC units]", 500, 300, 4000, 500, 300, 4000);

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
		h2_sci->Fill(sci21.E, sci22.E);
		h2_sci3v2->Fill(sci22.E, sci31.E);

		if(!mnd::IsInside(sci21.E, sci21_cut)) continue;
		if(!mnd::IsInside(sci22.E, sci22_cut)) continue;
		if(!mnd::IsInside(sci31.E, sci31_cut)) continue;

		h1_sci21_cut->Fill(sci21.E);
		h1_sci22_cut->Fill(sci22.E);
		h1_sci31_cut->Fill(sci31.E);

		for(const auto& cl : foot->fCl) {
			double delta = cl.Delta(); 
			double e = cl.fCE;
			double x = cl.fCX;
			
			/* Take the gain matching param. */
			e /= p->gain.CorrectionFactor(x);

			if(!mnd::IsInside(e, foot_cut)) continue;
			
			sum_energy_vs_delta -> Fill(delta, e);
			h2_m_vs_delta -> Fill(delta, cl.fCM);
			h1_delta -> Fill(delta);
			if(cl.fCM == 1) {
				h1_m1->Fill(e);
			} else if(cl.fCM > 2) {
				h1_m2->Fill(e);
			}
			h1_m0->Fill(e);
		}
	}

	/* Perform the Step (1) of the delta correction.
	 * Of course, this part isn't threadsafe. */
	auto* d = GetDeltaParams(*sum_energy_vs_delta);

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];
		if(sci21.hits.size() != 1) continue;
		if(sci22.hits.size() != 1) continue;
		if(sci31.hits.size() != 1) continue;

		if(sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut)) continue;
		if(sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut)) continue;
        if(sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut)) continue;

		for(const auto& cl : foot->fCl) {
			double delta = cl.Delta(); 
			double e = cl.fCE;
			double x = cl.fCX;

			/* Normalize again with the gain matching param. */
			e /= p->gain.CorrectionFactor(x);

			/* Normalize with the triangular delta-correction */
			e /= d->CorrectionBasic(delta);
			
			corr_energy_vs_delta -> Fill(delta, e / FOOTGainParam::CARBON_ADC - 1);
		}
	}
	
	/* Step (2) of the delta correction. */
	const double gaus_side_ratio = 1.4;
	Verbosity v = Verbosity::INFO;

#ifdef USING_FFT
	ROOT::DisableImplicitMT();	
	const int Q = 6;
	
	/* If: f(x) ≈ a0 + sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * Then e2 := e1 - sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * to fit everything onto the DC line. */
	auto [fft, gr, g] = DoFFTW<fit_info::GAUSS_MAX>(*corr_energy_vs_delta, -0.5, 0.5, gaus_side_ratio, Q, v);
	WARN("\n[[FOOT%d]] FFT Params:\n", p->de10_index_);
	printf("\"s\": %.5f,\n\"f\": %.5f,\n", d->s, d->f);
	printf("\"f2\": {\n\t\"n\": %d,", D_BINS);
	printf("\n\t\"c\": [");
	for(int i=0; i<Q; ++i) std::cout << fft.coeff[i] << ", ";
	std::cout << fft.coeff[Q] << "]";
	printf("\n}\n");

	ROOT::EnableImplicitMT();
#else
	constexpr size_t SPLINE_RANK = 12;
	WARN("Doing the polynomial fit with rank: %zu\n", SPLINE_RANK);
	auto [r, gr, g] = FitSplineAndGraph<SPLINE_RANK, fit_info::GAUSS_MAX>(*corr_energy_vs_delta, -0.5, 0.5, 80, gaus_side_ratio, v);
	std::cout << "Fit result e2(delta)  = " << r << std::endl;
#endif

	TCanvas *c = new TCanvas(Form("cRAW%d", ifoot), Form("Delta%d", ifoot), 2000, 1400);
	c->Divide(3,2);
	c->cd(1);
	h1_delta->Draw();
	c->cd(2);
	h2_m_vs_delta->Draw("COLZ");
	c->cd(3); gPad->SetLogz();
	sum_energy_vs_delta->Draw("COLZ");
	c->cd(4);  gPad->SetLogz();
	corr_energy_vs_delta->Draw("COLZ");
	g->Draw("L SAME");
	gr->Draw("P SAME");
	c->cd(5);
	h1_m1->Draw();
	c->cd(6);
	h1_m2->Draw();

	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* cs2 = new TCanvas("cs2", "SCI-Corr", 1600, 1200);
	cs2->Divide(2,1);
	cs2->cd(1); h2_sci->Draw("COLZ");
	cs2->cd(2); h2_sci3v2->Draw("COLZ");
}

FOOTDeltaParam* GetDeltaParams(TH2D* e_vs_delta) {
	constexpr static std::array<double, 2> SIDE_PEAK = {0.2, 0.4};
	double sp, sm, A0, Ap, Am;
	int bmax_0, bmax_p, bmax_m;
	auto* result = new FOOTDeltaParam{};

	auto h1_delta = std::unique_ptr<TH1D>( e_vs_delta->ProjectionX("px_") ); 
	h1_delta->SetDirectory(nullptr);
	
	TAxis* ax = h1_delta->GetXaxis();

	ax->SetRangeUser( -SIDE_PEAK[0], SIDE_PEAK[0] );
	bmax_0 = h1_delta->GetMaximumBin();
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( SIDE_PEAK[0], SIDE_PEAK[1] );
	bmax_p = h1_delta->GetMaximumBin();
	sp = ax->GetBinCenter(bmax_p);
	ax->SetRange(0,0); // unzoom

	ax->SetRangeUser( -SIDE_PEAK[1], -SIDE_PEAK[0] );
	bmax_m = h1_delta->GetMaximumBin();
	sm = ax->GetBinCenter(bmax_m);
	ax->SetRange(0,0); // unzoom
	double s = (sp - sm)/2;
	result->s = s;

	/* Next, in the TH2D just try to project the fews bins around sp,sm. */

	auto h1_e_0 = std::unique_ptr<TH1D>( e_vs_delta->ProjectionY("py_0", bmax_0-1, bmax_0+1) );
	auto h1_e_p = std::unique_ptr<TH1D>( e_vs_delta->ProjectionY("py_p", bmax_p-1, bmax_p+1) );
	auto h1_e_m = std::unique_ptr<TH1D>( e_vs_delta->ProjectionY("py_m", bmax_m-1, bmax_m+1) );
	h1_e_0->SetDirectory(nullptr);
	h1_e_p->SetDirectory(nullptr);
	h1_e_m->SetDirectory(nullptr);

	double side_ratio = 1;
	A0 = GaussFitMax( h1_e_0.get(), side_ratio, Verbosity::CHATTY ).first.at(1);
	Ap = GaussFitMax( h1_e_p.get(), side_ratio, Verbosity::CHATTY ).first.at(1);
	Am = GaussFitMax( h1_e_m.get(), side_ratio, Verbosity::CHATTY ).first.at(1);
	/* A0, by definition must be FOOTGainParam::CARBON_ADC, usually is within 1% */

	double f = (Ap + Am) / (2 * A0);
	result->f = f;

	WARN("\nA0 = %.2f, A+ = %.2f, A- = %.2f\n"
	       "s0 = %.2f, s+ = %.2f, s- = %.2f\n"
		   ">> f = %.4f; s = %.4f <<\n"
		   "b0 = %d, b+ = %d, b- = %d\n",
		   A0,Ap,Am, 0.0,sp,sm, f,s, bmax_0,bmax_p,bmax_m);

	return result;
}
