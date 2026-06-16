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

constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;

constexpr static std::array<double, 3> dq_binning = {120, -2.1, 2.1};
void foot_track_analysis (
	std::string fileName = "",
	std::vector <
		std::pair<u32, std::array<double, 2>>
	> cut_q = {},
	std::array<double, 3> binning_x = {100,-5,5},
	std::array<double, 3> binning_y = {100,-5,5},
	std::array<double, 3> binning_q = {120, 0.0, 7.5},
	DoSave do_save = DoSave::no
) {
	ROOT::EnableImplicitMT();

	TClass* cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!cl or !cl->GetDataMember("_x")) 
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	auto model = RNTupleModel::Create();
	auto foot = model->MakeField<RNFOOTHit>("FOOT");
	auto ntuple = RNTupleReader::Open(std::move(model), "h104", fileName);
	
	ROOT::EnableImplicitMT();

	TH1P *resx[N_PAIRS], *resy[N_PAIRS], *resq[N_PAIRS], *h1q[N_PAIRS], *h1_footx[N_PAIRS], *h1_footy[N_PAIRS];
	for(u32 i=0; i<N_PAIRS; ++i) {
		resx[i] = new TH1P(Form("((h1_fitrx_%d))Fit residue [mm]@FOOT%d X", i, i),
			kMagenta-9, binning_x[0], binning_x[1], binning_x[2]);
		resy[i] = new TH1P(Form("((h1_fitry_%d))Fit residue [mm]@FOOT%d Y", i, i),
			kYellow-9, binning_y[0], binning_y[1], binning_y[2]);
		resq[i] = new TH1P(Form("((h1_fitrq_%d))Fit residue [charge unit]@FOOT%d Q", i, i),
			kCyan-9, dq_binning[0], dq_binning[1], dq_binning[2]);

		h1q[i] = new TH1P(Form("((h1q_%d))Q value@FOOT%d", i, i),
			kBlue-9, binning_q[0], binning_q[1], binning_q[2]);
		h1_footx[i] = new TH1P(Form("((h1_footx_%d))FOOT%d X [mm]", i, i),
			kMagenta-9, 320, -50, 50);
		h1_footy[i] = new TH1P(Form("((h1_footy_%d))FOOT%d Y [mm]", i, i),
			kMagenta-9, 320, -50, 50);
	}
	TH2P* h2_resx = new TH2P("((h2_fitrx))Fit residue:FOOT ID@X-orientation",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_x[0], binning_x[1], binning_x[2]);
	TH2P* h2_resy = new TH2P("((h2_fitry))Fit residue:FOOT ID@Y-orientation",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_y[0], binning_y[1], binning_y[2]);
	TH2P* h2_resq = new TH2P("((h2_fitrq))Fit residue:FOOT ID@Y-orientation",
		N_PAIRS, -0.5, N_PAIRS-0.5, dq_binning[0], dq_binning[1], dq_binning[2]);

	TH2P* h2q = new TH2P("((h2q))Q value:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, binning_q[0], binning_q[1], binning_q[2]);

	TH2P* h2_kq = new TH2P("Kq cost [charge^2]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 500, 0, 1);
	TH2P* h2_kr = new TH2P("Kr cost [mm^2]:FOOT ID",
		N_PAIRS, -0.5, N_PAIRS-0.5, 400, 0, 50);
	for(const auto& cut : cut_q) {
		if(cut.first >= N_PAIRS)
			ERROR("Supplied index: %u >= %u as Q-cut pair index.\n", cut.first, N_PAIRS);
	}

	std::array<double, 3> xf, yf, zf, qf;

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);

		for(const RNFOOTTrack& t : foot->track) {
			if(t.n != N_PAIRS) continue;

			/* Check if veto is passed. */
			bool is_valid = true;
			for(const auto& cut : cut_q) {
				double qv = t._q[ cut.first ];
				if(!mnd::IsInside(qv, cut.second)) is_valid = false;
			}
			if(!is_valid) continue;

			for(size_t i=0; i<N_PAIRS; ++i) {
				double z = t._z[i];
				
				int w=0;
				/* Populate the arrays for N-1 fit. */
				for(size_t j=0; j<N_PAIRS; ++j) {
					if(j == i) continue;
					xf[w] = t._x[j];
					yf[w] = t._y[j];
					zf[w] = t._z[j];
					qf[w] = t._q[j];
					++w;
				}
				if(w != 3) ERROR("w must be 3, no?");

				auto fx = PolyFit<1>(zf, xf);
				auto fy = PolyFit<1>(zf, yf);
				//auto fx = std::array<double, 2> { t.x0, t.ax };	
				//auto fy = std::array<double, 2> { t.y0, t.ay };

				double x_extr = fx[0] + fx[1]*z;
				double y_extr = fy[0] + fy[1]*z;
				auto [q_extr, track_qvar] = mnd::mean_var(qf);
				
				double dx = t._x[i] - x_extr;
				double dy = t._y[i] - y_extr;
				double dq = t._q[i] - q_extr;
				double sq = t._sq[i]; // is sigma == sqrt(var);

				h1_footx[i]->Fill( t._x[i] );
				h1_footy[i]->Fill( t._y[i] );
				resx[i]->Fill(dx);
				resy[i]->Fill(dy);
				resq[i]->Fill(dq);
				h2_resx->Fill(i, dx);
				h2_resy->Fill(i, dy);
				h2_resq->Fill(i, dq);

				h1q[i]->Fill(t._q[i]);
				h2q->Fill(i, t._q[i]);
				
				/* Find the 3-measurement 'test' track charge params.. */
				h2_kq->Fill( i, sqrt(dq*dq + sq*sq) );
				h2_kr->Fill( i, sqrt(dx*dx + dy*dy) );
			}
		}
	}

	/* Around each residue, also fit a teeny-weeny gauss-chan 🥺 👉👈 */
	TCanvas* c1d = new TCanvas("FitResidue1D", "Fit residues 1D", 2150, 1650);
	c1d->Divide(N_PAIRS, 5);
	for(size_t i=0; i<N_PAIRS; ++i) {
		c1d->cd(i+1);
		auto [fit_result_x, _x_] = resx[i]->DrawAndFit(2.0, kMagenta + 2, 4.8);
		WARN("FOOT%zu X: " KBH_MAG "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_x[1], 1e3*fit_result_x[2]);
		c1d->cd(i+1 + N_PAIRS);
		auto [fit_result_y, _y_] = resy[i]->DrawAndFit(2.0, kMagenta + 2, 4.8);
		WARN("FOOT%zu Y: " KBH_YEL "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_y[1], 1e3*fit_result_y[2]);

		c1d->cd(i+1 + 2*N_PAIRS);
		auto [fit_result_q, _q_] = resq[i]->DrawAndFit(2.0, kMagenta + 2, 4.8);
		WARN("FOOT%zu Q: " KBH_CYN "%.3f, %.3f" KRNM " [Q]\n", i, fit_result_q[1], fit_result_q[2]);

		c1d->cd(i+1+ 3*N_PAIRS);
		h1_footx[i]->Draw();
		c1d->cd(i+1+ 4*N_PAIRS);
		h1_footy[i]->Draw();
	}
	
	TCanvas* c2d = new TCanvas("FitResidue2D", "Fit residues 2D", 2200, 1400);
	c2d->Divide(1, 3);
	c2d->cd(1); h2_resx->Draw("COLZ");
	c2d->cd(2); h2_resy->Draw("COLZ");
	c2d->cd(3); h2_resq->Draw("COLZ");

	TCanvas* cq1d = new TCanvas("Charge1D", "Charge of track measured by layers", 2000, 1200);
	cq1d->Divide(1, N_PAIRS);
	for(size_t i=0; i<N_PAIRS; ++i) {
		cq1d->cd(i+1);
		h1q[i]->Draw();
	}
	TCanvas* cq2d = new TCanvas("Charge2D", "Charge of track measured by layers", 2000, 1200);
	h2q->Draw("COLZ");
	
	TCanvas* ckv = new TCanvas("CostValues2D", "Cost Values 2D", 2200, 1400);
	ckv->Divide(2, 1);
	ckv->cd(1); h2_kq->Draw("COLZ");
	ckv->cd(2); h2_kr->Draw("COLZ");

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		save_all(canvas::Extension::png, { inf.stem().c_str() });
	}
}
