#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/PolyFitter.hxx"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/GaussFitMax.hxx"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

void tpc_alignment (
	std::string fileName = "", 
	int i_tpc=0,
	A3 binning_x = {100, -30, 30},
	A3 binning_y = {100, -30, 30},
	A2 sci21_cut = {-DBL_MAX, DBL_MAX},
	A2 sci22_cut = {-DBL_MAX, DBL_MAX},
	A2 sci31_cut = {-DBL_MAX, DBL_MAX},
	DoSave do_save = DoSave::no
) {
	if(i_tpc < 0 or i_tpc >= RNFRSCal::N_VALID_TPC)
		ERROR("Bad `i_tpc` arg. Must be [0,..%d]\n", RNFRSCal::N_VALID_TPC);

	const auto& label = RNFRSCal::tpc_label;
	WARN(BOLD ">>> Doing alignment on: TPC%s\n" KNRM, label[i_tpc]);
	ROOT::EnableImplicitMT();

	using Measurement = RNTPCCal::Measurement;
	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_params = f->Get < 
			std::remove_reference_t<decltype(*tpc_params)>
		> ("FRS_tpc_parameters");
		if(!tpc_params) ERROR("TPC param is nullptr.\n");
	}
	auto& tpc_param = tpc_params->at(i_tpc);

	TH1P *h1_x[2], *h1_y[4];
	for(int i=0; i<2; ++i)
		h1_x[i] = new TH1P(Form("TPC%d - delay line %d X [mm]@TPC%s,X", i_tpc, i, label[i_tpc]), kGreen-3, 
			binning_x[0], binning_x[1], binning_x[2]);

	for(int i=0; i<4; ++i)
		h1_y[i] = new TH1P(Form("TPC%d - anode %d Y [mm]@TPC%s,Y", i_tpc, i, label[i_tpc]), kRed-4,
			binning_y[0], binning_y[1], binning_y[2]);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

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

		auto& tpc = frs->tpc.at(i_tpc);

		const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;
		for(int d : {0,1}) {
			const std::vector<Measurement>& m = tpc_hits[d];
			if(m.size() != 1)
				continue;

			if(!std::isnan(m[0].x))
				h1_x[d]->Fill(m[0].x);

			for(int a : {0,1}) {
				if(!std::isnan(m[0].y[a]))
					h1_y[2*d + a]->Fill( m[0].y[a] );
			}
		}
	}

	std::cout << setprecision(10) << KRED;
	std::cout << "\"x_offset\": " << tpc_param.x_offset << std::endl;
	std::cout << "\"y_offset\": " << tpc_param.y_offset << KNRM << setprecision(6) << std::endl;
	/* Fit Gauss around those values. */
	for(int i=0; i<2; ++i) {
		auto [r, err] = GaussFitMax( *h1_x[i], 1.0); 
		auto [_, mu, sigma] = r;
		std::cout << "Fit result: \"DL" << i << ": " << r << "\n";
		tpc_param.x_offset[i] -= mu;
	}

	for(int i=0; i<4; ++i) {
		// Due to the TPC's y cutting out, no clue how wide is the gauss-chan 
		auto [r, err] = GaussFitMax( *h1_y[i], 1.0);
		auto [_, mu, sigma] = r;
		std::cout << "Fit result: \"AN" << i << ": " << r << "\n";
		tpc_param.y_offset[i] -= mu;
	}
	std::cout << setprecision(10) << KBH_GRN;
	std::cout << "\"x_offset\": " << tpc_param.x_offset << ",\n";
	std::cout << "\"y_offset\": " << tpc_param.y_offset << KNRM << ",\n";

	TCanvas* c = new TCanvas("TPC", "TPC", 2050, 1400);
	c->Divide(3,2);
	for(int d : {0,1}) {
		c->cd(3*d+1);
		h1_x[ d ]->Draw();
		for(int a: {0,1}) {
			c->cd(3*d + 2 + a);
			h1_y[ 2*d + a ]->Draw();
		}
	}
	
	TCanvas* cs = new TCanvas("SCIe", "SCI21,22,31", 1800, 800);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	
	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str(), Form("TPC%s", label[i_tpc]) });
	}
}
