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
#include "../../includes/util/MacroHelpers.hxx"

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

enum class Take { gauss, profile, gauss_fit_only };
enum class ShowOld { no, yes };

constexpr int D_BINS = 1000;
constexpr int N_STRIPS = 64;
constexpr double LR_SCALING = 2.5;

void foot_gain_match (
	std::string fileName = "", 
	int ifoot = 0,
	int bins_per_asic = 64,
	DoFit do_fit = DoFit::no,
	std::array<double,3> foot_binning = {1000, 4, 4000}, 
	double delta_cut = 0.05,
	uint32_t mult_cut = 1, // any multiplicity below that is disallowed 
	int Q_target = 6,
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no,
	Take take = Take::gauss,
	ShowOld show_old = ShowOld::no
) {
	if(N_STRIPS % bins_per_asic != 0)
		throw std::runtime_error(Form("Passed: %d , not evenly divisible by %d.\n",
			bins_per_asic, N_STRIPS));

	if(Q_target < 0 || Q_target > 6) 
		throw std::runtime_error(Form("Passed: %d, which isn't in {0,..=6}.\n", Q_target));

	const std::array<double, 2> delta_cut_mid   = { -delta_cut, delta_cut };

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
	auto* h1_delta_cut_mid = new TH1P("((h1_d1_mid))Delta [from -0.5, 0.5]", kRed-7, 150, -0.499, 0.499); 

	auto* h1_foot_e_mid = new TH1P("((h1_e_mid))FOOT E [ADC units]@Central strip value", ORGB{0xB2FD30}, (int)(1.5*foot_binning[0]), foot_binning[1], foot_binning[2]); 

	auto* hit_energy_mid = new TH2P(Form("((h2_mid))Cluster energy [ADC]:Strip number [0..640]@FOOT%d Raw, Requested Q=%d", ifoot, Q_target), 
		bins_per_asic*10, 0,640, foot_binning[0], foot_binning[1], foot_binning[2]);

	const double TARGET_ADC = Q_target*Q_target * 100.0;

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];
		
		if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
		if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
		if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);

		if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
		if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
		if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

		if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
		if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
		if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);

		for(const auto& cl : foot->fCl) {
			if(cl.fCM < mult_cut) continue;

			double delta = cl.Delta();
			h1_delta->Fill(delta);
			
			int i = static_cast<int>( cl.fCX );
			double e = cl.fCE;
				
			if(mnd::IsInside(delta, delta_cut_mid)) {
				h1_delta_cut_mid->Fill(delta);
				hit_energy_mid->FillInside(i, e);
				h1_foot_e_mid->FillInside(e); 
			}
		}
	}

	/* Idea is the following. Gain isn't always the same,.. some strips require higher gain for
	 * lower values. Simply to line up the total cluster energy values to the `TARGET_ADC`, across the detector. */
	
	/* Do the small fit in the 1D plot. */
	constexpr double sratio = 0.9;
	
	TH1D* h = *h1_foot_e_mid; 
	auto [fitr_mid, err_mid] = GaussFitMax(h, sratio);
	printf("[CENTRAL HIT] 1D projection yields: max: %.2f, gauss fit max (around this max+-%.1f sigma): %.2f +- %.2f\n",
		h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ), sratio,
		fitr_mid[1], err_mid[1]);
	
	const double mean_mid = fitr_mid[1];

	printf("\nCentral strip hit ADC gets mapped to: %s%.2f%s [ADC] since supplied: %sZ=%d%s\n", 
		BOLD, TARGET_ADC, KNRM, BOLD, Q_target, KNRM); 
	std::vector<TGraph*>      profile_fit, gauss_fit;
	std::vector<TGraphErrors*> profile_raw, gauss_raw;

	const std::vector<int> fit_these_asics = {2,3,4,5,6,7};
	auto contains = [](const auto& v, const typename std::decay_t<decltype(v)>::value_type& val) -> bool {
		return std::find(v.begin(), v.end(), val) != v.end();
	};

	constexpr static size_t POLY_DEG = 4;
	constexpr static int N_NEEDED_ENTRIES = 400;
	constexpr static int N_LOWEST_ENTRIES = 10;

	/* Try to fit a spline(s) for middle few ASICs. */
	if(do_fit == DoFit::yes) { 
		FOOTGainParam pp = foot_param->gain;

		const auto& v = do_fit.as_yes()->values;
		std::vector<double> fit_params;

		for(int a=0; a < TFOOTMapCont::N_ASIC; ++a) {
			FOOTAsicGainParam& asic = pp.fit.at(a);
			
			FMultiPoly* mp = asic.GetPoly(Q_target);
			if(!mp) {
				asic.multi_poly.emplace_back(Q_target); 
				mp = &asic.multi_poly.back();
			}
			
			double x_lo  = (a) * 64 + 0.00001;
			double x_hi = (a+1) * 64 - 0.00001;
			
			if( contains(v, a) ) {
				auto [rg, graw, gfit] = FitSplineAndGraph<POLY_DEG, fit_info::GAUSS_MAX> ( 
					*hit_energy_mid, x_lo, x_hi, 40, sratio /*, Verbosity::CHATTY */
				); 
				auto [rp, praw, pfit] = FitSplineAndGraph<POLY_DEG, fit_info::PROFILE_MAX> ( 
					*hit_energy_mid, x_lo, x_hi, 40 /*, 1.1, Verbosity::CHATTY */ 
				); 
				gauss_fit.push_back(gfit); 
				gauss_raw.push_back(graw); 
				profile_fit.push_back(pfit); 
				profile_raw.push_back(praw); 
				
				if(take == Take::gauss or take == Take::gauss_fit_only) {
					mp->pol = std::vector<double>(rg.begin(), rg.end()); 
				} else {
					mp->pol = std::vector<double>(rp.begin(), rp.end());
				}
			} 
			else {
				TAxis *xax = (*hit_energy_mid)->GetXaxis(); 
				int firstbin = xax->FindBin(x_lo); 
				int lastbin = xax->FindBin(x_hi); 

				auto pasic = std::unique_ptr<TH1D>((*hit_energy_mid)->ProjectionY("__py", firstbin, lastbin)); 
				pasic->SetDirectory(nullptr); 

				double profile_mean, gauss_mean; 
				if(pasic->Integral() >= N_NEEDED_ENTRIES) { /* If it contains more than 500 events, we can sample it. */ 
					profile_mean = pasic->GetXaxis()->GetBinCenter( pasic->GetMaximumBin() ); 
					auto [pg0, err_pg0] = GaussFitMax( pasic.get(), sratio );  
					gauss_mean = pg0[1]; 
				} else { /* No clue. Just take profile mean the mean, but gauss is invalidated. */ 
					profile_mean = pasic->GetMean(); 
					gauss_mean = mean_mid; 
				}
				TGraph* gfit = new TGraph(60); 
				TGraph* pfit = new TGraph(60); 
				for(int i=0; i<60; ++i) { 
					double x = x_lo + (i+0.00001) * (x_hi - x_lo)/59; 
					gfit->SetPoint(i, x, gauss_mean); 
					pfit->SetPoint(i, x, profile_mean); 
				} 
				if(pasic->Integral() >= N_NEEDED_ENTRIES) { 
					gfit->SetLineColor(gCol_); gfit->SetLineWidth(4); 
					pfit->SetLineColor(pCol_); pfit->SetLineWidth(4); 
				} else { 
					gfit->SetLineColor(gCol_ + 1); gfit->SetLineWidth(12); 
					pfit->SetLineColor(pCol_ + 1); pfit->SetLineWidth(12); 
				} 
				gauss_fit.push_back(gfit); 
				profile_fit.push_back(pfit); 

				mp->pol = std::vector<double>(1);
				if(pasic->Integral() < N_NEEDED_ENTRIES) {
					if(take == Take::gauss) { mp->pol[0] = gauss_mean; } 
					else if(pasic->Integral() >= N_LOWEST_ENTRIES) { mp->pol[0] = profile_mean; }
					else { mp->pol[0] = gauss_mean; }; // in case if there are really no entries there... 
				}
				else { // integral >= N_NEEDED_ENTRIES
					mp->pol[0] = gauss_mean; 
				}
			}

			std::sort( asic.multi_poly.begin(), asic.multi_poly.end() );
		}

		std::cout << "\"gain\": " << nlohmann::json(pp).dump(4) << std::endl;
		printf("Average value: %.5f (bin-center) and %.5f (gauss-fit-center)\n", 
			(*h1_foot_e_mid)->GetBinCenter((*h1_foot_e_mid)->GetMaximumBin()), mean_mid);
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
	
	TCanvas *c = new TCanvas("RawFOOT", Form("FOOT%d central", ifoot), 2400, 1400);
	hit_energy_mid->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1((*hit_energy_mid)->GetYaxis()->GetXmin());
		l->SetY2((*hit_energy_mid)->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	for(auto* pfit : profile_fit) pfit->Draw("L SAME");
	for(auto* praw : profile_raw) praw->Draw("P SAME");
	for(auto* gfit : gauss_fit) gfit->Draw("L SAME");
	for(auto* graw : gauss_raw) graw->Draw("P SAME");
	
	auto l = new TLegend(0.1,0.75,0.4,0.9);
	l->AddEntry(*hit_energy_mid, "Non-gain matched cluster energy (central hits)");
	l->AddEntry(gauss_fit[0], Form("Gaussian fit +-%.1f sigma around peak", sratio));
	l->AddEntry(profile_fit[0], "TProfile fit");
	if(show_old == ShowOld::yes) {
		const auto& gain = foot_param->gain;
		auto [g, _] = gain.GetRefZGraph(Q_target);
		g->SetLineColor(kPink - 2);
		g->Draw("L SAME");
		l->AddEntry(g, "Current gain curve from the setup file");
	}
	gStyle->SetLegendTextSize(0.021); l->Draw();

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	TCanvas* cd = new TCanvas("delta_energy", "Delta & 1D energy", 1400, 800);
	cd->Divide(2,2);
	cd->cd(1); h1_delta->Draw();
	cd->cd(2); h1_delta_cut_mid->Draw();
	cd->cd(3); h1_foot_e_mid->Draw();
	cd->cd(4);
	PLatex(0.08,
		Form("Requested charge: Q=%d", Q_target),
		Form("Detector ID: FOOT%d", ifoot),
		Form("Delta cut: #pm%.3f", delta_cut),
		Form("Cluster size >= %d", mult_cut),
		(take == Take::gauss) ? Form("Gauss fits size ratio: %.2f sigma", sratio)
			: "Fit taken from profile (violet curve)"
	);

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		if(do_fit == DoFit::yes) {
			save_all(canvas::Extension::png, 
				{
					Form("FOOT%d", ifoot),
					Form("Z_%d", Q_target),
					inf.stem().c_str()
				}
			);
		} else {
			save_all(canvas::Extension::png, 
				{
					Form("FOOT%d", ifoot),
					Form("Z_%d", Q_target),
					inf.stem().c_str(),
					"verification"
				}
			);
		}
	}
}
