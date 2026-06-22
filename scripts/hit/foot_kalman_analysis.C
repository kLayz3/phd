#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/Geometry.h"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

void foot_kalman_analysis (
	std::string fileName = "",
	DoSave do_save = DoSave::no
) {
	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTHit>("FOOT");
	auto frs = model->MakeField<RNFRSHit>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h104", fileName);
	
	constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	double Cr, Cq, Ct, Cp, max_cost;
	{
		std::array<double, 4>* c;
		TParameter<double>* m;
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2); Cp = c->at(3);
		max_cost = m->GetVal();
	}

	ROOT::EnableImplicitMT();

	
	TH1P* h1_track_mult = new TH1P("Track multiplicity", kRed-1, 10, -0.5, 9.5);
	TH2P* h2_q_vs_mult = new TH2P("Track charge [Q]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 40, 0,8);
	TH2P* h2_score_vs_mult = new TH2P("Track score [a.u.]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 300, 0, 50);
	
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		const size_t N = foot->track.size();
		h1_track_mult->Fill(N);

		for(const auto& t : foot->track) {
			h2_q_vs_mult->Fill(N, t.Q);
			h2_score_vs_mult->Fill(N, t.score);
		}
	}

	TCanvas* cm = new TCanvas("Multp", "Recognized tracks", 2150, 1400);
	cm->Divide(2,2);
	cm->cd(1); h2_q_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(3); h2_score_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(4); h1_track_mult->Draw();
	cm->cd(2); PLatex(0.08,
		"Coefficients: ",
		Form("Cr = %.1f mm^-2", Cr),
		Form("Cq = %.1f Q^-2", Cq),
		Form("Ct = %.1f mm^-2", Ct),
		Form("Cp = %.1f", Cp),
		Form("max cost [@4th]: %.1f", max_cost)
	);

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str() });
	}
}
