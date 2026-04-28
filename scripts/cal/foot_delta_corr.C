#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/GaussFitMax.hxx"
#include "../../includes/util/MacroHelpers.hxx"
#include "../../includes/util/FitSpline.hxx"

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
#endif

using namespace ROOT;
using namespace ROOT::Experimental;

struct DoFit {
	struct No {};
	struct Yes {
		std::vector<int> values;
	};

	/* constexpr */ static inline No no{};
	/* constexpr */ static inline Yes yes{};

	DoFit(No) : data_(No{}) {}
	DoFit(Yes y) : data_(std::move(y)) {}

	friend bool operator==(const DoFit& lhs, const DoFit& rhs) {
		return ( lhs.data_.index() == rhs.data_.index() &&
			lhs.data_.index() != std::variant_npos
		);
	}
	const Yes* as_yes() const { return std::get_if<Yes>(&data_); }
	
private:
	std::variant<No, Yes> data_;
};

/* Because cling issues the WEIRDEST compiler error,..
 * I cannot just simply return `FOOTDeltaParam` instance. It tries to compile the class
 * from the inputs but loses it on template spec in boost preproc library LOL. PEGI 18. */
class FOOTDeltaParam;
FOOTDeltaParam* GetDeltaParams(TH2D*);

constexpr double D_LO = 0.4999;
enum class PlotGain { no, yes };
enum class PlotRaw  { no, yes };

void foot_delta_corr (
	std::string fileName = "", 
	int ifoot = 0,
	std::array<double,3> foot_cut = {80, 1600,4600}, 
	int delta_bins = 80,
	int Q_target = 6,
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no,
	PlotGain plot_gain = PlotGain::no,
	PlotRaw plot_raw = PlotRaw::no,
	DoFit do_fit = DoFit::no
) {
#ifndef __USING_LUSTRE_HPC__
	gSystem->Load ("libfftw3.so");
#endif

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, 0, 
				i * 64, 1000);
		line->SetLineColor(kBlack);
		line->SetLineStyle(2);
		line->SetLineWidth(2);
		vlines.push_back( line );
	}
#define DRAW_VLINES(name) \
	{ \
		TH2D* h2_ = name; \
		for(auto* l0 : vlines) { \
			TLine* l = dynamic_cast<TLine*>(l0->Clone()); \
			l->SetY1(h2_->GetYaxis()->GetXmin()); \
			l->SetY2(h2_->GetYaxis()->GetXmax()); \
			l->Draw("SAME"); \
		} \
		gPad->SetGridx(false); gPad->Update(); \
	}

	FOOTParam* p;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		p = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!p)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
	}

	const double CA = FOOTGainParam::PROTON_ADC * Q_target * Q_target;
	const auto& gain = p->gain;

	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH1P* h1_delta = new TH1P(Form("((h1_foot%d))Delta@FOOT%d", ifoot, ifoot), kGreen+1, delta_bins, -D_LO, D_LO);
	TH2P* h2_m_vs_delta = new TH2P(Form("((h2_footm%d))Cluster size:delta@FOOT%d", ifoot, ifoot), delta_bins, -D_LO, D_LO, 10, 0.5, 10.5);
	TH2P* sum_energy_vs_x = new TH2P(Form("((h2_footraw%d))Cluster sum [ADC]:Strip num@FOOT%d gain matched", ifoot, ifoot), 
		160,0,640,
		foot_cut[0], foot_cut[1], foot_cut[2]);
	TH2P* sum_energy_vs_delta = new TH2P(Form("((h2_footraw%d))Cluster sum [ADC]:Delta@FOOT%d gain matched", ifoot, ifoot), 
		delta_bins,-D_LO, D_LO,
		foot_cut[0], foot_cut[1], foot_cut[2]);
	TH2P* corr_energy_vs_delta = new TH2P(Form("((h2_footcorr%d))E1/%.1f - 1:Delta@FOOT%d", ifoot, CA, ifoot), 
		delta_bins, -D_LO, D_LO,
		foot_cut[0], foot_cut[1]/CA - 0.9, foot_cut[2]/CA - 0.9);

	TH2P* h2_raw;
	if(plot_raw == PlotRaw::yes) 
		h2_raw = new TH2P(Form("((h2_raw))Raw Cluster ADC:Strip num@FOOT%d non-gain matched", ifoot),
		160,0,640, (int)(foot_cut[0] * 2), foot_cut[1] / 10, foot_cut[2] / 2.0);

	TH1P* h1_m0 = new TH1P(Form("((h_0)) FOOT%d dE [ADC units]@All cluster sizes", ifoot), kYellow - 7, foot_cut[0], foot_cut[1], foot_cut[2]);

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
			if(cl.fCM == 1) continue;

			double delta = cl.Delta(); 
			double e = cl.fCE;
			double x = cl.fCX;
			if(plot_raw == PlotRaw::yes) 
				h2_raw->Fill(x, e);
			e *= gain.CorrectionFactor(x, e);
			
			/* Explicitly skip here so we don't deal with over/underflow bins */
			if(!mnd::IsInside(e, foot_cut)) continue; 

			sum_energy_vs_delta -> Fill(delta, e);
			sum_energy_vs_x -> Fill(x, e);
			h2_m_vs_delta -> Fill(delta, cl.fCM);
			h1_delta -> Fill(delta);
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
			if(cl.fCM == 1) continue;

			double delta = cl.Delta(); 
			double e = cl.fCE;
			double x = cl.fCX;
			e *= gain.CorrectionFactor(x, e);

			/* Normalize with the triangular delta-correction */
			e /= d->CorrectionBasic(delta);
			
			corr_energy_vs_delta -> Fill(delta, e / CA - 1);
		}
	}
	
	/* Step (2) of the delta correction. */
	const double gaus_side_ratio = 1.4;
	Verbosity v = Verbosity::INFO;

#ifdef USING_FFT
	ROOT::DisableImplicitMT();	
	const int Q = 8;
	
	/* If: f(x) ≈ a0 + sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * Then e2 := e1 - sum_k [ a_k cos(2πk(x-xmin)/P) + b_k sin(2πk(x-xmin)/P) ]
	 * to fit everything onto the DC line. */
	auto [fft, gr, g] = DoFFTW<fit_info::GAUSS_MAX>(*corr_energy_vs_delta, -0.5, 0.5, gaus_side_ratio, Q, v);
	WARN("\n[[FOOT%d]] FFT Params:\n", p->de10_index_);
	printf("\"s\": %.5f,\n\"f\": %.5f,\n", d->s, d->f);
	printf("\"f2\": {\n\t\"n\": %d,", delta_bins);
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
	g->Draw("L SAME"); gPad->SetLogz();
	gr->Draw("P SAME");
	c->cd(5);
	sum_energy_vs_x->Draw("COLZ");
	c->cd(6);
	h1_m0->Draw();

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

	if(plot_gain == PlotGain::yes) {
		TCanvas* gn = new TCanvas("gain", "Gain factor", 2000, 1300);
		gn->Divide(2,2);
		TH2D* hgain = gain.GetHisto();
		hgain->SetTitle(Form("FOOT%d Gain Parameter", ifoot));
		gn->cd(1); hgain->Draw("COLZ"); DRAW_VLINES(hgain)
		gn->cd(2); hgain->Draw("SURF1");
		auto [p4, graph4] = gain.GetGraph(4*64 + 32);
		gn->cd(3); graph4->Draw("AL"); p4->Draw("P SAME");
		PLatex(0.08, "Middle of ASIC[4]"); 
		auto [p5, graph5] = gain.GetGraph(5*64 + 32);
		gn->cd(4); graph5->Draw("AL"); p5->Draw("P SAME");
		PLatex(0.08, "Middle of ASIC[5]"); 
	}
	if(plot_raw == PlotRaw::yes) {
		TCanvas* craw =  new TCanvas("non-gain-matched", "Non gain matched", 1400, 800);
		h2_raw->Draw("COLZ"); DRAW_VLINES(*h2_raw)
		TGraph* ref6 = gain.GetRefZGraph(6);
		TGraph* ref5 = gain.GetRefZGraph(5);
		TGraph* ref4 = gain.GetRefZGraph(4);
		TGraph* ref3 = gain.GetRefZGraph(3);
		ref5->SetLineColor(kMagenta + 1);
		ref4->SetLineColor(kPink - 2);
		ref3->SetLineColor(kOrange + 7);
		ref6->Draw("L SAME");
		ref5->Draw("L SAME");
		ref4->Draw("L SAME");
		ref3->Draw("L SAME");

		auto l = new TLegend(0.1,0.75,0.38,0.9);
		l->AddEntry(*h2_raw, "Non-gain matched cluster energy");
		l->AddEntry(ref6, "Z=6 reference");
		l->AddEntry(ref5, "Z=5 reference");
		l->AddEntry(ref4, "Z=4 reference");
		l->AddEntry(ref3, "Z=3 reference");
		gStyle->SetLegendTextSize(0.027);
		l->Draw();
		
		if(do_fit == DoFit::yes and Q_target == 6) {
			// Try to fit low energy (high-delta part)
			const auto& v = do_fit.as_yes()->values;

			auto contains = [](const auto& v, const typename std::decay_t<decltype(v)>::value_type& val) -> bool {
				return std::find(v.begin(), v.end(), val) != v.end();
			};
			constexpr static size_t POLY_DEG = 4;
			constexpr double sratio = 0.1;
			for(int a : v) {
				double x_lo  = (a) * 64 + 0.00001;
				double x_hi = (a+1) * 64 - 0.00001;
			
				auto [rg, graw, gfit] = FitSplineAndGraph<POLY_DEG, fit_info::GAUSS_MAX> ( 
					*h2_raw, x_lo, x_hi, 40, sratio /*, Verbosity::CHATTY */
				); 
				gfit->Draw("L SAME");
				graw->Draw("P SAME");
				printf("ASIC[%d]: ", a); std::cout << rg << std::endl;
				const auto* z4 = gain.fit[a].GetPoly(4);
				FMultiPoly z4_new {};
				z4_new.Z = 4;
				z4_new.pol = std::vector(rg.begin(), rg.end());

				std::cout << "Old: " << nlohmann::json(*z4).dump(4) << std::endl;
				std::cout << "New: " << nlohmann::json(z4_new).dump(4) << std::endl;
			}
		}
	}

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { Form("FOOT%d", ifoot), inf.stem().c_str() });
	}
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
