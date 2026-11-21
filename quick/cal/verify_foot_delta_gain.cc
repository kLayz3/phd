#include "util/CLI.h"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/TextUtils.h"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

#include "monad/accessories/formatting.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

constexpr f64 DELTA_EDGE_DEFAULT = 0.4999;

template<typename T>
struct Cont {
	T x,y;
	T& operator[](u32 n) noexcept {
		switch(n) {
			case(0): return x;
			case(1): return y;
			default: ERROR("Requested {} != 0,1 in {}\n", n, MND_SELF_TYPE_CSTR);
		}
	}
	const T& operator[](u32 n) const noexcept {
		switch(n) {
			case(0): return x;
			case(1): return y;
			default: ERROR("Requested {} != 0,1 in {}\n", n, MND_SELF_TYPE_CSTR);
		}
	}
};

using TH2Pair = Cont<TH2P>;
using TH1Pair = Cont<TH1P>;
using FPair   = Cont<std::shared_ptr<RNFOOTCal>>;

int main(int argc, char* argv[]) {
	CLI::App app{"Verify the gain-matching and δ-correction of a specific FOOT detector pair."};

	std::vector<std::string> fileNames {};
	i32 ipair = 0;
	A3 foot_x_raw_binning = {80, 1600,4600};
	A3 foot_y_raw_binning = {80, 1600,4600};
	A3 foot_x_gain_binning = {80, 1600,4600};
	A3 foot_y_gain_binning = {80, 1600,4600};
	A3 delta_binning = {80, -DELTA_EDGE_DEFAULT, DELTA_EDGE_DEFAULT};
	
	u32 nbins_pos = 640;
	
	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	A3 sci_binning = {500, 300, 4000};
	
	std::vector<canvas::Extension> save{};
	std::string desc;
	bool plot_gain = false;

	u32 nthreads = 1;

	add_logged_option(app, "-f,--file", fileNames, "Pass a file name.")
		->check(CLI::ReadPermissions)
        ->delimiter(',');
	add_logged_option(app, "-i,--foot-id", ipair,
		"Select which FOOT detector pair.")
		->check(CLI::Range(0, N_FOOT_DETECTORS/2 - 1));
	add_logged_option(app, "-n,--nbins", nbins_pos, "Select number of bins along the detector position axis.");
	add_logged_option(app, "--foot-x-raw-binning", foot_x_raw_binning, "FOOT Raw ADC binning (X-measuring detector)")
		->delimiter(',');
	add_logged_option(app, "--foot-y-raw-binning", foot_y_raw_binning, "FOOT Raw ADC binning (Y-measuring detector)")
		->delimiter(',');
	add_logged_option(app, "--foot-x-gain-binning", foot_x_gain_binning, "FOOT Gain ADC binning (X-measuring detector)")
		->delimiter(',');
	add_logged_option(app, "--foot-y-gain-binning", foot_y_gain_binning, "FOOT Gain ADC binning (Y-measuring detector)")
		->delimiter(',');
	add_logged_option(app, "-d,--delta-binning", delta_binning, "FOOT delta-1D binning")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "--sci-binning",sci_binning, "Scintillator QDC binning")
		->delimiter(',');
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');
	add_logged_option(app, "-desc", desc, "Folder name where to save. Ignored if there is no -o,--save passed.");
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
	)->check(CLI::PositiveNumber);
	add_logged_flag(app, "-g,--plot-gain,!--no-plot-gain", plot_gain, "Do plot gain");
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.size() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

	Cont<FOOTParam*> p;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileNames.front().c_str(), "READ");
		get_obj(f, p.x, Form("FOOT%d_setup", 2*ipair));
		get_obj(f, p.y, Form("FOOT%d_setup", 2*ipair + 1));
	}
	if(p.x->GetOrientation() == p.y->GetOrientation())
		ERROR("Both detectors in the pair pointing at identical orientation?");
	if(p.x->GetOrientation() == Orientation::UNKNOWN
	|| p.y->GetOrientation() == Orientation::UNKNOWN)
		ERROR("One of the detectors have unknown orientation?");

	const u32 ix = (p.x->GetOrientation() == Orientation::X) ? 2*ipair : 2*ipair+1;
	const u32 iy = (p.x->GetOrientation() == Orientation::Y) ? 2*ipair : 2*ipair+1;
	const A3 pos_binning = { (double)nbins_pos, 0, RNFOOTCal::N_STRIPS };

	ROOT::EnableThreadSafety();
	ROOT::DisableImplicitMT();
	TApplication rootApp("app", 0, 0);

	TH2Pair raw_energy_vs_x {
		.x = TH2P{Form("((h2_raw_x))Cluster E [ADC gain-adjusted]:x [strip #]@FOOT%u raw-spectrum (x-measuring)", ix),
			pos_binning, foot_x_raw_binning},
		.y = TH2P{Form("((h2_raw_y))Cluster E [ADC gain-adjusted]:y [strip #]@FOOT%u raw-spectrum (y-measuring)", iy),
			pos_binning, foot_y_raw_binning},
	};
	TH2Pair gain_energy_vs_delta {
		.x = TH2P{Form("((h2_gain_x))Cluster E [ADC gain-adjusted]:Delta@FOOT%u gain-matched (x-measuring)", ix),
			delta_binning, foot_x_gain_binning},
		.y = TH2P{Form("((h2_gain_y))Cluster E [ADC gain-adjusted]:Delta@FOOT%u gain-matched (y-measuring)", iy),
			delta_binning, foot_y_gain_binning},
	};
	TH2Pair gain_energy_vs_x {
		.x = TH2P{Form("((h2_gain_x))Cluster E [ADC gain-adjusted]:x [strip #]@FOOT%u gain-matched (x-measuring)", ix),
			pos_binning, foot_x_gain_binning},
		.y = TH2P{Form("((h2_gain_y))Cluster E [ADC gain-adjusted]:x [strip #]@FOOT%u gain-matched (y-measuring)", iy),
			pos_binning, foot_y_gain_binning
		},
	};
	TH1Pair gain_e {
		.x = TH1P{Form("((h1_gain_x))dE [ADC gain-adjusted]@FOOT%u gain-matched (x-measuring)", ix),
			kYellow - 5, foot_x_gain_binning},
		.y = TH1P{Form("((h1_gain_y))dE [ADC gain-adjusted]@FOOT%u gain-matched (y-measuring)", iy),
			kYellow - 5, foot_y_gain_binning},
	};
	TH2Pair final_energy_vs_delta {
		.x = TH2P{Form("((h2_final_x))Cluster E [ADC gain-delta adjusted]:Delta@FOOT%u gain-delta matched (x-measuring)", ix),
			delta_binning, foot_x_gain_binning},
		.y = TH2P{Form("((h2_final_y))Cluster E [ADC gain-delta adjusted]:Delta@FOOT%u gain-delta matched (y-measuring)", iy),
			delta_binning, foot_y_gain_binning},
	};
	TH2Pair final_energy_vs_x {
		.x = TH2P{Form("((h2_gain_x))Cluster E [ADC gain-delta adjusted]:x [strip #]@FOOT%u gain-delta matched (x-measuring)", ix),
			pos_binning, foot_x_gain_binning},
		.y = TH2P{Form("((h2_gain_y))Cluster E [ADC gain-delta adjusted]:x [strip #]@FOOT%u gain-delta matched (y-measuring)", iy),
			pos_binning, foot_y_gain_binning
		},
	};
	TH1Pair final_e {
		.x = TH1P{Form("((h1_gain_x))dE [ADC gain-delta adjusted]@FOOT%u gain-delta matched (x-measuring)", ix),
			kYellow - 7, foot_x_gain_binning},
		.y = TH1P{Form("((h1_gain_y))dE [ADC gain-delta adjusted]@FOOT%u gain-delta matched (y-measuring)", iy),
			kYellow - 7, foot_y_gain_binning},
	};
	TH1Pair cl_mult {
		.x = TH1P{Form("((h1_cl_mult_x))Cluster multiplicity@FOOT%u (x-measuring)", ix),
			kCyan -3, 10, 0.5, 10.5},
		.y = TH1P{Form("((h1_cl_mult_y))Cluster multiplicity@FOOT%u (y-measuring)", iy),
			kCyan -3, 10, 0.5, 10.5},
	};
	TH2Pair cl_mult_vs_delta {
		.x = TH2P{Form("((h2_cl_mult_x))Cluster multiplicity:Delta@FOOT%u (x-measuring)", ix),
			delta_binning, 10, 0.5, 10.5},
		.y = TH2P{Form("((h2_cl_mult_y))Cluster multiplicity:Delta@FOOT%u (y-measuring)", iy),
			delta_binning, 10, 0.5, 10.5},
	};

	auto h2_both = TH2P{Form("((h2_0)) FOOT%u dE x-measuring [ADC corrected]:FOOT%d dE y-measuring [ADC corrected]",
		ix, iy), foot_x_gain_binning, foot_y_gain_binning};
	
	auto h1_sci21 = TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, sci_binning);
	auto h1_sci22 = TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, sci_binning);
	auto h1_sci31 = TH1P("SCI31 QDC mean [QDC units]", 0x009B2F_c, sci_binning);
	auto h1_sci21_cut = TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, sci_binning);
	auto h1_sci22_cut = TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, sci_binning);
	auto h1_sci31_cut = TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, sci_binning);

	show_console_cursor(false);

	DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
		auto model = RNTupleModel::Create();
		FPair foot;
		foot.x = model->MakeField<RNFOOTCal>(Form("FOOT%u", ix));
		foot.y = model->MakeField<RNFOOTCal>(Form("FOOT%u", iy));
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fname);

        const size_t nentries = ntuple->GetNEntries();

		auto _bar_ptr = std::make_unique<ProgressBar>(
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{fmt::format(
					"{}/{}: {:L} ({})", i+1, fileNames.size(), nentries, fname
				)},
                option::ForegroundColor{ indicators::next_col() },
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];
		
		constexpr size_t MHIT_CAPACITY = 20;
		Cont<f64> delta, x, ei, ef;
		Cont<std::array<f64, MHIT_CAPACITY>> hits_e;
		Cont<size_t> nhits = {.x = 0, .y = 0};

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			
			const auto& sci21 = frs->sci[0];
			const auto& sci22 = frs->sci[1];
			const auto& sci31 = frs->sci[2];
			
			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() == 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() == 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() == 1) h1_sci31_cut->Fill(sci31.E);

#define EXPAND_(V) \
	delta.V  = cl##V.Delta(); \
	x.V  = cl##V.fCX; \
	ei.V = cl##V.fCE; \
	\
	ei.V *= p.V->gain.CorrectionFactor(x.V, ei.V); \
	\
	gain_energy_vs_delta.V -> Fill(delta.V, ei.V); \
	gain_energy_vs_x.V -> Fill(x.V, ei.V); \
	gain_e.V->Fill(ei.V); \
	\
	ef.V = ei.V / p.V->de.CorrectionFactor( delta.V ); \
	final_energy_vs_delta.V -> Fill(delta.V, ef.V); \
	final_energy_vs_x.V -> Fill(x.V, ef.V); \
	final_e.V->Fill( ef.V ); \
	\
	cl_mult.V -> Fill( cl##V.fCM ); \
	cl_mult_vs_delta.V -> Fill( delta.V, cl##V.fCM ); \
	/* Export the energy of this cluster to the outside world. */ \
	if(nhits.V < MHIT_CAPACITY) \
		hits_e.V[nhits.V++] = ef.V; \
		
			for(const auto& clx : foot.x->fCl) {
				EXPAND_(x)
			}
			for(const auto& cly : foot.y->fCl) {
				EXPAND_(y)
			}

			for(size_t nx = 0; nx < nhits.x; ++nx)
				for(size_t ny = 0; ny < nhits.y; ++ny)
					h2_both->Fill(hits_e.x[nx], hits_e.y[ny]);
		} // for(auto entryId : *ntuple)

		if(!bar.is_completed())
			bar.mark_as_completed();
    });

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
		TH2D* h2_ = &name.h; \
		for(auto* l0 : vlines) { \
			TLine* l = dynamic_cast<TLine*>(l0->Clone()); \
			l->SetY1(h2_->GetYaxis()->GetXmin()); \
			l->SetY2(h2_->GetYaxis()->GetXmax()); \
			l->Draw("SAME"); \
		} \
	}
	TCanvas* cg = new TCanvas("GainMatched", "Single detector (gain-matched)", 2100, 1400);
	TCanvas* cf = new TCanvas("FinalSpectra", "Single detector (gain-delta matched)", 2100, 1400);
	cg->Divide(3,2); cf->Divide(3,2);
	
	for(u32 i: {0,1}) {
		cg->cd(3*i + 1); gPad->SetLogz();
		gain_energy_vs_delta[i].Draw("COLZ");
		cg->cd(3*i + 2); gPad->SetLogz();
		gain_energy_vs_x[i].Draw("COLZ");
		DRAW_VLINES(gain_energy_vs_x[i]);

		cg->cd(3*i + 3); 
		gain_e[i].Draw();

		cf->cd(3*i + 1); gPad->SetLogz();
		final_energy_vs_delta[i].Draw("COLZ");
		cf->cd(3*i + 2); gPad->SetLogz();
		final_energy_vs_x[i].Draw("COLZ");
		DRAW_VLINES(final_energy_vs_x[i]);

		cf->cd(3*i + 3); 
		final_e[i].Draw();
	}

	TCanvas* c2 = new TCanvas("PairXY", Form("FOOT Pair %d XY", ipair), 2000, 1400);
	c2->Divide(2,1);
	c2->cd(1); gPad->SetLogz();
	h2_both.Draw("COLZ");
	c2->cd(2); new PLatex(0.08,
		Form("FOOT%u measuring '%s'", ix, p.x->orientation.c_str()),
		Form("FOOT%u measuring '%s'", iy, p.x->orientation.c_str()),
		Form("Cut applied on SCI21: (%.1f, %.1f)", sci21_cut[0], sci21_cut[1]),
		Form("Cut applied on SCI22: (%.1f, %.1f)", sci22_cut[0], sci22_cut[1])
	);
	
	TCanvas* cs = new TCanvas("SCI1D", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21.Draw();
	cs->cd(2); h1_sci22.Draw();
	cs->cd(3); h1_sci31.Draw();
	cs->cd(4); h1_sci21_cut.Draw();
	cs->cd(5); h1_sci22_cut.Draw();
	cs->cd(6); h1_sci31_cut.Draw();
	
	TCanvas* cmult = new TCanvas("multiplicities", "multiplicities", 2400, 1200);
	cmult->Divide(2,2);
	for(u32 i: {0,1}) {
		cmult->cd(2*i + 1);
		cl_mult[i].Draw();
		cmult->cd(2*i + 2);
		cl_mult_vs_delta[i].Draw("COLZ");
	}

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOTPair%d", ipair),
		desc,
		mnd::rootfile_names_concatenated(fileNames)
	});

	WARN("End-of-main");
	rootApp.Run();
	return 0;
}
