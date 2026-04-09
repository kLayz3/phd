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
#include "../../includes/util/Tracking.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

struct DoFit {
	struct No {};
	struct Verify {};
	struct Yes {
		std::vector<int> values;
	};

	/* constexpr */ static inline No no{};
	/* constexpr */ static inline Verify verify{};
	/* constexpr */ static inline Yes yes{};

	DoFit(No) : data_(No{}) {}
	DoFit(Verify) : data_(Verify{}) {}
	DoFit(Yes y) : data_(std::move(y)) {}

	friend bool operator==(const DoFit& lhs, const DoFit& rhs) {
		return ( lhs.data_.index() == rhs.data_.index() &&
			lhs.data_.index() != std::variant_npos
		);
	}
	const Yes* as_yes() const { return std::get_if<Yes>(&data_); }
	
private:
	std::variant<No, Verify, Yes> data_;
};

constexpr int D_BINS = 1000;
constexpr int N_STRIPS = 64;
constexpr double LR_SCALING = 2.5;

enum class HitPos { mid, lr, none };

void foot_gain_match (
	std::string fileName = "", 
	int ifoot = 0,
	int bins_per_asic = 64,
	DoFit do_fit = DoFit::no,
	std::array<double,3> foot_binning = {1000, 4, 4000}, 
	std::array<double,4> delta_cut = {0.27, 0.05, 0.03, 0.03},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX} 
) {
	if(N_STRIPS % bins_per_asic != 0)
		throw std::runtime_error(Form("Passed: %d , not evenly divisible by %d.\n",
			bins_per_asic, N_STRIPS));
	
	const std::array<double, 2> delta_cut_mid   = { -delta_cut[3], delta_cut[3] };
	const std::array<double, 2> delta_cut_left = { -delta_cut[0]-delta_cut[1], -delta_cut[0]+delta_cut[2] };
	const std::array<double, 2> delta_cut_right  = { delta_cut[0]-delta_cut[2], delta_cut[0]+delta_cut[1] };

	FOOTParam *foot_param; 
	TParameter<bool> *is_already_gain_matched;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		foot_param = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!foot_param)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
		is_already_gain_matched = f->Get<TParameter<bool>>(Form("FOOT%d_gain_matched", ifoot));
		if(!is_already_gain_matched)
			ERROR("FOOT%d gain matching tag not fetchable? Line: %d\n", ifoot, __LINE__); 
	}
	if(do_fit == DoFit::yes and is_already_gain_matched->GetVal()) 
		ERROR("Trying to fit the gain matched curve, but FOOT%d tagged as already gain matched!\n"
			"Either use this to verify, or re-do the file without gain-matching flag upfront.", ifoot);
	
	if(do_fit == DoFit::verify and !is_already_gain_matched->GetVal()) 
		ERROR("Trying to verify the gain matched data, but FOOT%d tagged as NOT already gain matched!\n"
			"Either use this to gain-match, or re-do the file WITH gain-matching flag upfront.", ifoot);

	ROOT::EnableImplicitMT();
	
	auto get_hit_position = [&delta_cut_mid, &delta_cut_right, &delta_cut_left](double d) { 
		if(mnd::IsInside(d, delta_cut_mid)) return HitPos::mid;
		if(mnd::IsInside(d, delta_cut_left) || mnd::IsInside(d, delta_cut_right)) return HitPos::lr;
		return HitPos::none;
	};

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
	auto* h1_delta_cut_mid = new TH1P("((h1_d1_mid))Delta [from -0.5, 0.5]", kRed-7, 150, -0.499, 0.499); 
	auto* h1_delta_cut_lr = new TH1P("((h1_d1_lr))Delta [from -0.5, 0.5]", kBlue-7, 150, -0.499, 0.499); 

	auto* h1_foot_e_mid = new TH1P("((h1_e_mid))FOOT E [ADC units]@Central strip value", ORGB{0xB2FD30}, (int)(1.5*foot_binning[0]), foot_binning[1], foot_binning[2]); 
	auto* h1_foot_e_lr = new TH1P("((h1_e_lr))FOOT E [ADC units]@Central strip value", ORGB{0x4AFFFF}, (int)(1.5*foot_binning[0]), foot_binning[1]/LR_SCALING, foot_binning[2]/LR_SCALING); 

	auto* hit_energy_mid = new TH2P(Form("((h2_mid))Max Signal [ADC]:Strip number [0..640]@FOOT%d Raw, per ASIC", ifoot), 
		bins_per_asic*10, 0,640, foot_binning[0], foot_binning[1], foot_binning[2]);
	auto* hit_energy_lr = new TH2P(Form("((h2_lr))Max Signal [ADC]:Strip number [0..640]@FOOT%d Raw, per ASIC", ifoot), 
		bins_per_asic*10, 0,640, foot_binning[0], foot_binning[1]/LR_SCALING, foot_binning[2]/LR_SCALING);

	const double C_ADC = FOOTGainParam::CARBON_ADC;
	std::array<double, 3> v_binning = {foot_binning[0], C_ADC/(3*LR_SCALING), C_ADC*1.5};

	auto* hit_energy_corr = new TH2P(Form("((h2_foot%d_corr))Max strip ADC Corrected [A.U.]:Strip number [0..640]@FOOT%d Corrected, per ASIC", ifoot, ifoot), 
		bins_per_asic*10, 0,640, v_binning[0], v_binning[1], v_binning[2]);
	auto* h1_hit_energy_corr = new TH1P(Form("((h1_foot%d_corr))Max strip ADC Corrected [A.U.]@FOOT%d Corrected, per ASIC", ifoot, ifoot), ORGB{0x52FD30}, v_binning[0], v_binning[1], v_binning[2]); 
	auto* clust_energy_corr = new TH2P(Form("((h2_foot%d_corr))Cluster Sum Corrected [ADC]:Strip number [0..640]@FOOT%d Corrected, per ASIC", ifoot, ifoot), 
		bins_per_asic*10, 0,640, v_binning[0], v_binning[1], v_binning[2]);
	auto* h1_clust_energy_corr = new TH1P(Form("((h1_foot%d_corr))Cluster Sum Corrected [A.U.]@FOOT%d Corrected, per ASIC", ifoot, ifoot), ORGB{0x52FD30}, v_binning[0], v_binning[1], v_binning[2]); 

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
			
			int i = static_cast<int>( cl.fCX );
			double e = cl.fCP; /* Take *central* strip value, not the entire cluster! */
		
			auto hp = get_hit_position(delta);
			switch(hp) {
				case HitPos::mid: {
					h1_delta_cut_mid->Fill(delta);
					hit_energy_mid->Fill(i, e);
					h1_foot_e_mid->Fill(e); break;
				}
				case HitPos::lr: {
					h1_delta_cut_lr->Fill(delta);
					hit_energy_lr->Fill(i, e);
					h1_foot_e_lr->Fill(e); break;
				}
				default: break;
			}
			if(do_fit == DoFit::verify /* && hp != HitPos::none */) {
				double cog = cl.fCX;
				double e_cl = cl.fCE;
				hit_energy_corr->Fill(cog, e);
				h1_hit_energy_corr->Fill(e);

				clust_energy_corr->Fill(cog, e_cl);
				h1_clust_energy_corr->Fill(e_cl);
			}
		}
	}

	/* Idea is the following. Gain isn't always the same,.. some strips require higher gain for
	 * lower values. Simply to line up both the lateral delta and central ones, across the detector.
	 * Simple single factor-per-strip gain matching won't cut it here. */
	
	/* Do the small fit in the 1D plot. */
	const double sratio = 1.1;
	TH1D* h = *h1_foot_e_mid; 
	auto [fitr_mid, err_mid] = GaussFitMax(h, sratio);
	printf("[CENTRAL] 1D projection yields: max: %.2f, gauss fit max (around this max+-%.1f sigma): %.2f +- %.2f\n",
		h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ), sratio,
		fitr_mid[1], err_mid[1]);
	h = *h1_foot_e_lr; 
	auto [fitr_lr, err_lr] = GaussFitMax(h, sratio);
	printf("[LATERAL] 1D projection yields: max: %.2f, gauss fit max (around this max+-%.1f sigma): %.2f +- %.2f\n",
		h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ), sratio,
		fitr_lr[1], err_lr[1]);
	
	double S,M;
	double mean_mid = M = fitr_mid[1];
	double mean_lr  = S = fitr_lr[1];
	double C_ADC_LR = C_ADC * S / M;

	printf("\nIf central gets mapped to: %s%.2f%s [ADC] then lateral is: %s%.2f%s\n\n", BOLD, C_ADC, KNRM, BOLD, C_ADC_LR, KRNM); 
	std::vector<TGraph*>      profile_fit_mid, gauss_fit_mid;
	std::vector<TGraphErrors*> profile_raw_mid, gauss_raw_mid;
	std::vector<TGraph*>      profile_fit_lr, gauss_fit_lr;
	std::vector<TGraphErrors*> profile_raw_lr, gauss_raw_lr;

	const std::vector<int> fit_these_asics = {2,3,4,5,6,7};
	auto contains = [](const auto& v, const typename std::decay_t<decltype(v)>::value_type& val) -> bool {
		return std::find(v.begin(), v.end(), val) != v.end();
	};

	constexpr static size_t POLY_DEG = 4;
	constexpr static int N_NEEDED_ENTRIES = 250;

	/* I'm coding like absolute monkey. Legit tilted fuckin hell... */
#define EXPAND_FIT(pre) \
	{ \
		auto [rg, graw, gfit] = FitSplineAndGraph<POLY_DEG, fit_info::GAUSS_MAX> ( \
			*hit_energy_##pre, x_lo, x_hi, 40, 1.1 /*, Verbosity::CHATTY */ \
		); \
		auto [rp, praw, pfit] = FitSplineAndGraph<POLY_DEG, fit_info::PROFILE_MAX> ( \
			*hit_energy_##pre, x_lo, x_hi, 40 /*, 1.1, Verbosity::CHATTY */ \
		); \
		gauss_fit_##pre.push_back(gfit); \
		gauss_raw_##pre.push_back(graw); \
		profile_fit_##pre.push_back(pfit); \
		profile_raw_##pre.push_back(praw); \
\
		fit_params_##pre = std::vector<double>(rg.begin(), rg.end()); \
	}

#define EXPAND_PROJ(pre) \
	{ \
		TAxis *xax = (*hit_energy_##pre)->GetXaxis(); \
		int firstbin = xax->FindBin(x_lo); \
		int lastbin = xax->FindBin(x_hi); \
 \
		auto pasic = std::unique_ptr<TH1D>((*hit_energy_##pre)->ProjectionY("__py", firstbin, lastbin)); \
		pasic->SetDirectory(nullptr); \
 \
		double profile_mean, gauss_mean; \
		if(pasic->Integral() > N_NEEDED_ENTRIES) { /* If it contains less than 500 events, just skip it. */ \
			profile_mean = pasic->GetXaxis()->GetBinCenter( pasic->GetMaximumBin() ); \
			auto [pg0, err_pg0] = GaussFitMax( pasic.get(), 1.0 );  \
			gauss_mean = pg0[1]; \
		} else { /* No clue. Just take profile mean the same, but gauss is invalidated. */ \
			profile_mean = pasic->GetXaxis()->GetBinCenter( pasic->GetMaximumBin() ); \
			gauss_mean = mean_##pre; \
		} \
		TGraph* gfit = new TGraph(60); \
		TGraph* pfit = new TGraph(60); \
		for(int i=0; i<60; ++i) { \
			double x = x_lo + (i+0.00001) * (x_hi - x_lo)/59; \
			gfit->SetPoint(i, x, gauss_mean); \
			pfit->SetPoint(i, x, profile_mean); \
		} \
		if(pasic->Integral() > N_NEEDED_ENTRIES) { \
			gfit->SetLineColor(gCol_); gfit->SetLineWidth(4); \
			pfit->SetLineColor(pCol_); pfit->SetLineWidth(4); \
		} else { \
			gfit->SetLineColor(gCol_ + 1); gfit->SetLineWidth(12); \
			pfit->SetLineColor(pCol_ + 1); pfit->SetLineWidth(12); \
		} \
		gauss_fit_##pre.push_back(gfit); \
		profile_fit_##pre.push_back(pfit); \
\
		fit_params_##pre = std::vector<double>(1); \
		fit_params_##pre[0] = gauss_mean; \
	}
	/* Try to fit a spline(s) for middle few ASICs. */
	
	if(do_fit == DoFit::yes) { 
		auto h2_gain = new TH2P(Form("Gain factor:Raw FOOTE [ADC]@FOOT%d for different strips", ifoot),
			2000, 0, foot_binning[2], 1000, 0.0, 5.0);
		auto h2_gain_relative = new TH2P(Form("Gain factor relative:Raw FOOTE [ADC]@FOOT%d for different strips", ifoot),
			2000, 0, foot_binning[2], 1000, 0.5, 3.0);

		TGraph* cutoff = new TGraph(); 
		TGraph* hi_pts = new TGraph(); TGraph* hi_pts1 = new TGraph();
		TGraph* lo_pts = new TGraph(); TGraph* lo_pts1 = new TGraph();

		FOOTGainParam pp {};
		pp.mid_avg = M;
		pp.lat_avg = S;

		//const auto& v = fit_these_asics;
		const auto& v = do_fit.as_yes()->values;
		std::vector<double> fit_params_mid, fit_params_lr;

		for(int a=0; a < TFOOTMapCont::N_ASIC; ++a) {
			double x_lo  = (a) * 64 + 0.00001;
			double x_hi = (a+1) * 64 - 0.00001;
			if( contains(v, a) ) {
				EXPAND_FIT(mid)
				EXPAND_FIT(lr)
			}
			else {
				EXPAND_PROJ(mid)
				EXPAND_PROJ(lr)
			}
			pp.fit[a].central = fit_params_mid;
			pp.fit[a].lateral = fit_params_lr;
			//std::cout << "[MID]: FOOT" <<  foot_param->de10_index_ << ": parameters (ASIC: " << a << "): " << fit_params_mid << std::endl;
			//std::cout << "[LR ]: FOOT" <<  foot_param->de10_index_ << ": parameters (ASIC: " << a << "): " << fit_params_lr << std::endl;

			static constexpr int N_PTS_G = 64;
			double dx = (x_hi - x_lo) / (N_PTS_G - 1);
			for(int i = 0; i<N_PTS_G; ++i) {
				double x = x_lo + (i+0.5) * dx;
				double phi0 = poly::Eval(x, fit_params_mid);
				double phi1 = poly::Eval(x, fit_params_lr);
				
				double gain0 = C_ADC / phi0;
				double gain1 = C_ADC_LR / phi1;
				auto gainLine = GetLine( {phi0, gain0}, {phi1, gain1} );
				auto gainLineRelative = GetLine (  {phi0, 1}, {phi1, gain1/gain0} );
				const auto& [offset, slope] = gainLine;
				
				double gain_times_e_derivative_zero = - (offset / (2*slope));
				if( std::isfinite(gain_times_e_derivative_zero) ) {
					cutoff->AddPoint( gain_times_e_derivative_zero, offset / 2);
				}
				hi_pts ->AddPoint(phi0, gain0);
				lo_pts ->AddPoint(phi1, gain1);
				hi_pts1->AddPoint(phi0, 1);
				lo_pts1->AddPoint(phi1, gain1/gain0);

				FillTrack(*h2_gain, gainLine);
				FillTrack(*h2_gain_relative, gainLineRelative);
			}
		}
		cutoff->SetMarkerStyle(20);
		cutoff->SetMarkerSize(1.1);
		cutoff->SetMarkerColor(kMagenta +1);

		hi_pts ->SetMarkerStyle(22);
		hi_pts ->SetMarkerSize(1.15);
		hi_pts ->SetMarkerColor(kYellow -3);
		hi_pts1->SetMarkerStyle(22);
		hi_pts1->SetMarkerSize(1.15);
		hi_pts1->SetMarkerColor(kYellow -3);

		lo_pts ->SetMarkerStyle(23);
		lo_pts ->SetMarkerSize(1.15);
		lo_pts ->SetMarkerColor(kOrange +4);
		lo_pts1->SetMarkerStyle(23);
		lo_pts1->SetMarkerSize(1.15);
		lo_pts1->SetMarkerColor(kOrange +4);
		
		TCanvas *cgain = new TCanvas("cgain", "Gain Factor", 2200, 1400);
		cgain->Divide(1,2);

		cgain->cd(1); h2_gain->Draw("COLZ"); (*h2_gain)->SetStats(0); 
		cutoff->Draw("P SAME");
		hi_pts ->Draw("P SAME");
		lo_pts ->Draw("P SAME");
		TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
		leg->AddEntry(hi_pts, "Central strip hit region gain", "p");
		leg->AddEntry(lo_pts, "Side strip hit region gain", "p");
		leg->AddEntry(cutoff, "Cutoff : #frac{d(g*e)}{de} = 0", "p");
		leg->Draw();

		cgain->cd(2); (*h2_gain_relative)->SetStats(0); h2_gain_relative->Draw("COLZ");
		hi_pts1->Draw("P SAME");
		lo_pts1->Draw("P SAME");
		TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
		leg1->AddEntry(hi_pts1, "Central strip hit region gain", "p");
		leg1->AddEntry(lo_pts1, "Side strip hit region gain", "p");
		leg1->Draw();

		std::cout << "\"gain\": " << nlohmann::json(pp).dump(4) << std::endl;
	}

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, foot_binning[1], 
				                i * 64, foot_binning[2]);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}
	
	TCanvas *c = new TCanvas(Form("cRAW%d", ifoot), Form("FOOT%d central", ifoot), 2000, 1400);
	hit_energy_mid->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1((*hit_energy_mid)->GetYaxis()->GetXmin());
		l->SetY2((*hit_energy_mid)->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	
	for(auto* pfit : profile_fit_mid) pfit->Draw("L SAME");
	for(auto* praw : profile_raw_mid) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit_mid) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw_mid) graw->Draw("P SAME");

	TCanvas *clr = new TCanvas(Form("cRAW_lr%d", ifoot), Form("FOOT%d peripheral", ifoot), 2000, 1400);
	hit_energy_lr->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1((*hit_energy_lr)->GetYaxis()->GetXmin());
		l->SetY2((*hit_energy_lr)->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	for(auto* pfit : profile_fit_lr) pfit->Draw("L SAME");
	for(auto* praw : profile_raw_lr) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit_lr) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw_lr) graw->Draw("P SAME");

	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	auto* hs_d = new THStack("hs_d","Delta cut");
	hs_d->Add(*h1_delta_cut_mid);
	hs_d->Add(*h1_delta_cut_lr);
	TCanvas* cd = new TCanvas("cd", "Delta & 1D energy", 1400, 800);
	cd->Divide(2,2);
	TLine* l_dl = new TLine(-delta_cut[0], (*h1_delta)->GetXaxis()->GetXmin(), 
							-delta_cut[0], (*h1_delta)->GetXaxis()->GetXmin());
	TLine* l_dr = new TLine( delta_cut[0], (*h1_delta)->GetXaxis()->GetXmin(), 
							 delta_cut[0], (*h1_delta)->GetXaxis()->GetXmin());
	l_dl->SetLineColor(kRed); l_dr->SetLineColor(kRed);
	l_dl->SetLineStyle(2);    l_dr->SetLineStyle(2);
	l_dl->SetLineWidth(3);    l_dr->SetLineWidth(3);
	cd->cd(1); h1_delta->Draw(); l_dl->Draw("SAME"); l_dr->Draw("SAME");
	cd->cd(2); hs_d->Draw();
	cd->cd(3); h1_foot_e_mid->Draw();
	cd->cd(4); h1_foot_e_lr->Draw();

	if(do_fit == DoFit::verify) {
		TCanvas *cc = new TCanvas(Form("cCORR%d", ifoot), Form("Corrected FOOT%d", ifoot), 2000, 1400);
		cc->Divide(2,2);
#define DRAW_AND_DO_LINES \
		h2->Draw("COLZ"); \
		for(auto* l0 : vlines) { \
			TLine* l = dynamic_cast<TLine*>(l0->Clone()); \
			l->SetY1(h2->GetYaxis()->GetXmin()); \
			l->SetY2(h2->GetYaxis()->GetXmax()); \
			l->Draw("SAME"); \
		} 
		TH2D* h2 = &hit_energy_corr->h;
		cc->cd(1);
		DRAW_AND_DO_LINES

		h2 = &clust_energy_corr->h;
		cc->cd(2);
		DRAW_AND_DO_LINES

		cc->cd(3);
		h1_hit_energy_corr->Draw();
		cc->cd(4);
		h1_clust_energy_corr->Draw();
	}
}
