#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

constexpr static double NBINS  = 300;
constexpr static double CUT_LO = 10;
constexpr static double CUT_HI = 3000;

enum class DoMagic {no, yes};

void verify_foot_delta_gain (
	std::string fileName = "", 
	int i_pair = 0,
	std::array<double,3> cut0 = {NBINS, CUT_LO, CUT_HI}, 
	std::array<double,3> cut1 = {NBINS, CUT_LO, CUT_HI},
	std::array<double,2> sci21_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci22_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> sci31_cut = {-DBL_MAX, DBL_MAX},
	std::array<double,2> delta_axis = {80, 0.501},
	uint32_t mult_cut = 1, // every cluster with size below this gets rejected
	DoMagic do_magic = DoMagic::no,
	DoSave do_save = DoSave::no
) {
	if(i_pair > 3) throw std::runtime_error("i_pair parameter (2nd argument) is > 3.");
	ROOT::EnableImplicitMT();

	std::array<FOOTParam*, 2> p;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		p[0] = f->Get<FOOTParam>(Form("FOOT%d_setup", 2*i_pair));
		p[1] = f->Get<FOOTParam>(Form("FOOT%d_setup", 2*i_pair + 1));
		if(!p[0] or !p[1])
			throw std::runtime_error(Form("FOOT param is nullptr. Fix it (line: %d).", __LINE__));
	}

	std::array<std::shared_ptr<RNFOOTCal>, 2> foot {};
	auto model = RNTupleModel::Create();
	for(int i: {0,1}) 
		foot[i] = model->MakeField<RNFOOTCal>(Form("FOOT%d", 2*i_pair + i));
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH2P* gain_energy_vs_delta[2];
	TH2P* gain_energy_vs_x[2];
	TH1P* gain_e[2];
	TH2P* final_energy_vs_delta[2];
	TH2P* final_energy_vs_x[2];
	TH1P* final_e[2];
	TH1P* cl_mult[2];
	TH2P* cl_mult_vs_delta[2];
	TH2P* cl_mult_e[2][5];

	for(int i: {0,1}) {
		const auto& lim = (i==0) ? cut0 : cut1;
		int foot_i = 2*i_pair + i;
		
		/* ================================================================== */
		gain_energy_vs_delta[i] = new TH2P(Form("((h2_gain_%d))Cluster E [ADC gain-adjusted]:Delta@FOOT%d gain-matched", i, foot_i),
			delta_axis[0], -1.0*delta_axis[1], delta_axis[1], lim[0], lim[1], lim[2]	
		);
		gain_energy_vs_x[i] = new TH2P(Form("((h2_gain_x%d))Cluster E [ADC gain-adjusted]:x [strip #]@FOOT%d gain-matched", i, foot_i),
			640, 0, 640, lim[0], lim[1], lim[2]
		);
		gain_e[i] = new TH1P(Form("((h1_gain_e%d))dE [ADC gain-adjusted]@FOOT%d gain-matched", i, foot_i), 
			kYellow - 5, lim[0], lim[1], lim[2]
		);
		/* ================================================================== */
		final_energy_vs_delta[i] = new TH2P(Form("((h2_final_d%d))Cluster E [ADC gain-delta adjusted]:Delta@FOOT%d gain-delta adjusted", i, foot_i),
			delta_axis[0], -1.0*delta_axis[1], delta_axis[1], lim[0], lim[1], lim[2]	
		);
		final_energy_vs_x[i] = new TH2P(Form("((h2_final_x%d))Cluster E [ADC gain-delta adjusted]:x [strip #]@FOOT%d gain-delta adjusted", i, foot_i),
			640, 0, 640, lim[0], lim[1], lim[2]	
		);
		final_e[i] = new TH1P(Form("((h1_final_e%d))dE [ADC gain-delta adjusted]@FOOT%d gain-delta adjusted", i, foot_i), 
			kYellow - 7, lim[0], lim[1], lim[2]
		);
		cl_mult[i] = new TH1P(Form("((h1_cl_mult%d))Cluster multiplicity@FOOT%d", i, foot_i),
			kCyan -3, 10, 0.5, 10.5);
		cl_mult_vs_delta[i] = new TH2P(Form("((h2_cl_mult%d))Cluster multiplicity:Delta@FOOT%d", i, foot_i),
			delta_axis[0], -1.0*delta_axis[1], delta_axis[1], 10, 0.5, 10.5);
		for(int m: {1,2,3,4,5}) {
			int index = m-1;
			cl_mult_e[i][index] = new TH2P(Form("((h2_cl_mult_e%d_%d))Cluster E [ADC gain-delta adjusted]:x [strip #]@FOOT%d, mult=%d", i, index, foot_i, m),
				640, 0, 640, lim[0], lim[1], lim[2]);
		}
	}

	TH2P* h2_both = new TH2P(Form("((h2_0)) FOOT%d dE [ADC corrected]:FOOT%d dE [ADC corrected]", 2*i_pair+1, 2*i_pair), 
		cut0[0], cut0[1], cut0[2],
		cut1[0], cut1[1], cut1[2]);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);
	auto* h2_sci2v2 = new TH2P("SCI22 QDC mean [QDC units]:SCI21 QDC mean [QDC units]", 500, 300, 4000, 500, 300, 4000);
	auto* h2_sci3v2 = new TH2P("SCI31 QDC mean [QDC units]:SCI22 QDC mean [QDC units]", 500, 300, 4000, 500, 300, 4000);
	auto* h2_sci2v2_cut = new TH2P("((h2_sc_cut))SCI22 QDC mean [QDC units]:SCI21 QDC mean [QDC units]@With cut", 500, 300, 4000, 500, 300, 4000);
	auto* h2_sci3v2_cut = new TH2P("((h2_sc_cut))SCI31 QDC mean [QDC units]:SCI22 QDC mean [QDC units]@With cut", 500, 300, 4000, 500, 300, 4000);

	double x[2]; // strip x
	double delta[2]; // delta
	double gfactor[2]; // gain factor
	double dfactor[2]; // delta factor
	double ei[2]; // ADC intermediate (gain-corrected)
	double ef[2]; // ADC final (gain-corrected + delta-corrected)

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		const auto& sci21 = frs->sci[0];
		const auto& sci22 = frs->sci[1];
		const auto& sci31 = frs->sci[2];

		if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
		if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
		if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);
		if(sci31.hits.size() >= 1 and sci22.hits.size() >= 1) 
			h2_sci3v2->Fill(sci22.E, sci31.E);
		if(sci21.hits.size() >= 1 and sci22.hits.size() >= 1) 
			h2_sci2v2->Fill(sci21.E, sci22.E);

		if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
		if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
		if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

		if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
		if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
		if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);
		if(sci31.hits.size() >= 1 and sci22.hits.size() >= 1) 
			h2_sci3v2_cut->Fill(sci22.E, sci31.E);
		if(sci21.hits.size() >= 1 and sci22.hits.size() >= 1) 
			h2_sci2v2_cut->Fill(sci21.E, sci22.E);

#define EXPAND_(i) \
	delta[i]  = cl##i.Delta(); \
	x[i]  = cl##i.fCX; \
	ei[i] = cl##i.fCE; \
	\
	ei[i] *= p[i]->gain.CorrectionFactor(x[i], ei[i]); \
	\
	gain_energy_vs_delta[i] -> Fill(delta[i], ei[i]); \
	gain_energy_vs_x[i] -> Fill(x[i], ei[i]); \
	gain_e[i]->Fill(ei[i]); \
	\
	dfactor[i] = p[i]->de.CorrectionFactor( delta[i] ); \
	ef[i] = ei[i] / dfactor[i]; \
	final_energy_vs_delta[i] -> Fill(delta[i], ef[i]); \
	final_energy_vs_x[i] -> Fill(x[i], ef[i]); \
	final_e[i]->Fill( ef[i] ); \
	\
	cl_mult[i] -> Fill( cl##i.fCM ); \
	cl_mult_vs_delta[i] -> Fill( delta[i], cl##i.fCM ); \
	switch(cl##i.fCM) { \
		case(1): { \
			cl_mult_e[i][0] -> Fill(x[i], ef[i]); break; \
		} \
		case(2): { \
			cl_mult_e[i][1] -> Fill(x[i], ef[i]); break; \
		} \
		case(3): { \
			cl_mult_e[i][2] -> Fill(x[i], ef[i]); break; \
		} \
		case(4): { \
			cl_mult_e[i][3] -> Fill(x[i], ef[i]); break; \
		} \
		case(5): { \
			cl_mult_e[i][4] -> Fill(x[i], ef[i]); break; \
		} \
		default: { break; } \
	}
			
		for(const auto& cl0 : foot[0]->fCl) {
			if(cl0.fCM < mult_cut) continue;
			EXPAND_(0)
			
			for(const auto& cl1 : foot[1]->fCl) {
				if(cl1.fCM < mult_cut) continue;
				EXPAND_(1)
				h2_both->Fill(ef[0], ef[1]);
			}
		}
	}
	
	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * 64, 0, 
				                i * 64, 1000);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}
#define DRAW_VLINES(name) \
	{ \
		TH2D* h2_ = &name->h; \
		for(auto* l0 : vlines) { \
			TLine* l = dynamic_cast<TLine*>(l0->Clone()); \
			l->SetY1(h2_->GetYaxis()->GetXmin()); \
			l->SetY2(h2_->GetYaxis()->GetXmax()); \
			l->Draw("SAME"); \
		} \
	}
	TCanvas* cg = new TCanvas("GainMatched", "Single (gain-matched)", 2100, 1400);
	TCanvas* cf = new TCanvas("FinalSpectra", "Single (final)", 2100, 1400);
	cg->Divide(3,2); cf->Divide(3,2);
	
	for(int i: {0,1}) {
		cg->cd(3*i + 1); gPad->SetLogz();
		gain_energy_vs_delta[i]->Draw("COLZ");
		cg->cd(3*i + 2); gPad->SetLogz();
		gain_energy_vs_x[i]->Draw("COLZ");
		DRAW_VLINES(gain_energy_vs_x[i]);

		cg->cd(3*i + 3); 
		gain_e[i]->Draw();

		cf->cd(3*i + 1); gPad->SetLogz();
		final_energy_vs_delta[i]->Draw("COLZ");
		cf->cd(3*i + 2); gPad->SetLogz();
		final_energy_vs_x[i]->Draw("COLZ");
		DRAW_VLINES(final_energy_vs_x[i]);

		cf->cd(3*i + 3); 
		final_e[i]->Draw();
	}

	TCanvas* c2 = new TCanvas("PairXY", Form("FOOT Pair %d XY", i_pair), 2000, 1400);
	c2->Divide(2,1);
	c2->cd(1); gPad->SetLogz();
	h2_both->Draw("COLZ");
	c2->cd(2);
	PLatex(0.08, 
		Form("FOOT%d measuring '%s'", 2*i_pair+0, p[0]->orientation.c_str()),
		Form("FOOT%d measuring '%s'", 2*i_pair+1, p[1]->orientation.c_str()),
		Form("Cut applied on SCI21: (%.1f, %.1f)", sci21_cut[0], sci21_cut[1]),
		Form("Cut applied on SCI22: (%.1f, %.1f)", sci22_cut[0], sci22_cut[1])
	);
	
	TCanvas* cs = new TCanvas("SCI1D", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	TCanvas* cs2 = new TCanvas("SCI2D", "SCI-Corr", 2200, 1400);
	cs2->Divide(2,2);
	cs2->cd(1); h2_sci2v2->Draw("COLZ");
	cs2->cd(2); h2_sci3v2->Draw("COLZ");
	cs2->cd(3); h2_sci2v2_cut->Draw("COLZ");
	cs2->cd(4); h2_sci3v2_cut->Draw("COLZ");


	TCanvas* cmult = new TCanvas("multiplicities", "multiplicities", 2400, 1200);
	cmult->Divide(5,2);
	for(int i: {0,1}) {
		cmult->cd(5*i + 1);
		cl_mult[i]->Draw();
		cmult->cd(5*i + 2);
		cl_mult_vs_delta[i]->Draw("COLZ");
		
		for(int m: {1,2,3}) {
			cmult->cd(5*i + 2 + m);
			cl_mult_e[i][m-1] -> Draw("COLZ");	
		}
	}

#define PARSE_CUT_ARGS(x) \
	(sci##x##_cut[0] > 0   ) ? std::to_string((int)sci##x##_cut[0]).c_str() : "inf", \
	(sci##x##_cut[1] < 5000) ? std::to_string((int)sci##x##_cut[1]).c_str() : "inf" 

	if(do_save == DoSave::yes) {
		std::filesystem::path inf( fileName );
		std::string postfix = std::string( Form("21_%s-%s;22_%s-%s;33-%s-%s",
			PARSE_CUT_ARGS(21), PARSE_CUT_ARGS(22), PARSE_CUT_ARGS(31)) );
		save_all(canvas::Extension::png, { Form("pair%d", i_pair), inf.stem().c_str(), postfix });
	}
}
