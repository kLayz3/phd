/* In this script, we align XY of each of the FOOT detectors
 * Alignment to FRS is done later by offseting the entire box. */

#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/FitSpline.hxx"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;
enum class DoVerify { no, yes };

/* Here we place FOOT ifoot onto various positions z, oriented either x- or y- 
 * and see which one fits the picture as correlation. */
void foot_spread (
	std::string fileName = "",
	uint32_t ifoot = 0, 
	std::array<double,3> binning_x = {200,-30,30},
	std::array<double,2> foot_q_cut = {5.4, 6.6},
	std::vector<std::pair<double, double>> dead_regions = {}, 
	std::array<double,2> sci21_cut = {NAN, NAN},
	std::array<double,2> sci22_cut = {NAN, NAN},
	std::array<double,2> sci31_cut = {NAN, NAN},
	DoSave do_save = DoSave::no,
	DoVerify do_verify = DoVerify::no
) {
	if(ifoot > 7 or ifoot == 4 or ifoot == 5)
		throw std::invalid_argument("Second argument `ifoot` can be only {0,1,2,3, 6,7}.");

	if(do_verify == DoVerify::yes)
		do_save = DoSave::no;

	FOOTParam *foot_param; 
	FOOTBoxParam *box;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		foot_param = f->Get<FOOTParam>(Form("FOOT%d_setup", ifoot));
		if(!foot_param)
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
		box = f->Get<FOOTBoxParam>("FOOT0_box");
		if(!box)
			throw std::runtime_error(Form("FOOT box param is nullptr. Fix it (line: %d).", __LINE__));
	}
	Orientation o = foot_param->GetOrientation();
	if(o == Orientation::UNKNOWN)
		throw std::runtime_error("FOOT orientation not specified. I won't allow it.\n");

	auto* foot_q_vs_d = new TH2P(Form("Cluster Charge:Delta [-0.5, 0.5]@FOOT%d", ifoot), 
		80, -0.5, 0.5, 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_q_vs_x = new TH2P(Form("Cluster Charge:FOOT measurement [mm]@FOOT%d", ifoot), 
		binning_x[0], binning_x[1], binning_x[2], 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_pos = new TH1P(Form("((h1))FOOT measurement [mm]@FOOT%d", ifoot), kCyan - 6, binning_x[0], binning_x[1], binning_x[2]);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);

	auto model = RNTupleModel::Create();
	auto frs  = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
	ROOT::EnableImplicitMT();

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

		/* In an event, only a single valid FOOT cluster must be found. */
		bool is_foot_event_valid = false;
		double xFOOT_ = NAN;
		double qFOOT_ = NAN;
		double dFOOT_ = NAN;
		for(const auto& hit : foot->fCl) {
			if(hit.fCM == 1) continue;	
			double q = foot_param->Q( hit );
			if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;
			
			double d = hit.Delta();
			double cx = hit.fCX;
			
			double hit_position;
			if(do_verify == DoVerify::yes)
				hit_position = foot_param->X0(hit);
			else
				hit_position = foot_param->BarePosition(hit);

			if( std::isfinite(xFOOT_) ) {
				/* Already found valid point in the event. */
				is_foot_event_valid = false; break;
			} else {
				/* Export it outside. */
				xFOOT_ = hit_position;
				qFOOT_ = q;
				dFOOT_ = d;
			}
			is_foot_event_valid = true;
		}
		if(!is_foot_event_valid) continue;

		foot_pos->Fill(xFOOT_);
		foot_q_vs_d->Fill(dFOOT_, qFOOT_);
		foot_q_vs_x->Fill(xFOOT_, qFOOT_);
	}
	/* Some strips/regions are dead and not to confuse the gauss-chan 🥺 👉👈,
	 * exclude these regions from the fit. Passed in as a region sequence. */
	constexpr double fit_area_factor = 0.5;
	auto mygaus = [dead_regions](Double_t* x, Double_t* par) -> Double_t {
		double xx = x[0];
		for(const auto& [lo, hi] : dead_regions) {
			if(xx >= lo && xx <= hi) {
				TF1::RejectPoint();
				return 0.0;
			}
		}
		return par[0] * std::exp(
			-0.5 * std::pow((xx-par[1])/par[2], 2)
		);
	};
	ROOT::DisableImplicitMT();

	TF1* f = new TF1("f", mygaus, fit_area_factor*binning_x[1], fit_area_factor*binning_x[2], 3);
	f->SetParNames("A", "mu", "sigma");
    f->SetParameters((*foot_pos)->GetMaximum(), (*foot_pos)->GetMean(), (*foot_pos)->GetStdDev());
	f->SetBit(TObject::kCanDelete, false);
	f->SetBit(TF1::kNotGlobal);

	(*foot_pos)->Fit(f, "RMN");
	WARN("Result to be put in setup file: " BOLD "\"delta_p\": %.5f\n" KNRM, -f->GetParameter("mu"));
	TCanvas* c = new TCanvas("FOOT", Form("FOOT%d calibration position plots", ifoot), 2150,1400);
	c->Divide(2,2);
	c->cd(1); foot_q_vs_d->Draw("COLZ"); gPad->SetLogz();
	c->cd(2); foot_q_vs_x->Draw("COLZ"); gPad->SetLogz();
	c->cd(3); foot_pos->Draw("COLZ");
	
	f->Draw("same");

	TCanvas* cs = new TCanvas("SCIs&TPCs", "SCI21,22,31 and TPC ref", 2200, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("FOOT%d", ifoot) });
	}
}
