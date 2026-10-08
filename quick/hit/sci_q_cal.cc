#include "TFRSCalCont.h"
#include "util/CLI.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"

#include "TApplication.h"
#include "TFRSHitCont.h"
#include "TFOOTHitCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

constexpr auto SCI21_I = RNFRSCal::SCI21_I;
constexpr auto SCI22_I = RNFRSCal::SCI22_I;
constexpr auto SCI31_I = RNFRSCal::SCI31_I;
constexpr auto SCI41_I = RNFRSCal::SCI41_I;
constexpr auto N_VALID_SCI = RNFRSCal::N_VALID_SCI;

int main(int argc, char* argv[]) {
	CLI::App app{"Calibrate the QDC into-charge measurement of SCI 21/22/31/41. "
		"Prerequisite before running this program is to have both the pedestal subtraction and "
        "velocity dependence (β) correction done. These programs are from the cal directory; "
		"`sci_qdc_ped` and `tof_cal`, respectively. "
		"Optionally, can also use FOOT's tracking here to imply the charge Q."};

	std::vector<std::string> fileName;
	u32 i_sci = 0;
	u32 nthreads = 1;
	mnd::Option<u32> q = mnd::None;
	double d = 0.4;
	A3 bin = {400,100,4000};
    u32 niter = 2;
    unsigned short line_size = 4;
    double sratio = 1.4;

	double beta_range = 0.2;
	constexpr int beta_binning = 200;

	add_logged_option(app, "-f,--file", fileName, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-i,--sci", i_sci,
        mnd::msg("Scintillator index; %u => SCI21, %u => SCI22, %u => SCI31, %u => SCI41",
            SCI21_I, SCI22_I, SCI31_I, SCI41_I))
        ->check(CLI::IsMember(
			{SCI21_I, SCI22_I, SCI31_I, SCI41_I}
		));
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-q,--charge", q, "Select the charge: 1...Qm on which to cut on FOOT's. "
		"If left as undefined, or in None state, won't cut on FOOT.");
	add_logged_option(app, "-d,--range", d, "Select the charge interval width, [q-d, q+d], on which to cut on FOOT's")
		->check(CLI::Range(0.01,0.50));
	add_logged_option(app, "--binning", bin, "Binning for SCI QDC axis.")
		->delimiter(',');
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
        ->check(CLI::PositiveNumber);
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
        ->check(CLI::PositiveNumber);
	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
	add_logged_option(app, "-b,--beta-range", beta_range, "β-interval around nominal β taken from Bρ")
		->check(CLI::Range(0.05, 0.5));

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;
	if(fileName.size() == 0)
		ERROR("To continue, must supply at least one file name!\n");
	
    const auto& label = RNFRSCal::sci_label;
	const mnd::Option<A2> qcut = q.map([d](u32 val) {
		return A2{val-d, val+d};
	});

	FRSToFParam const* tof_param {};
	BeamInfo const* beam_param {};
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params {};
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.front().c_str(), "READ");
        get_obj(f, sci_params, "FRS_sci_parameters");
        get_obj(f, tof_param,  "FRS_tof_par");
		get_obj(f, beam_param, "FRS_beam_info");
	}
	const SCIParam& par = sci_params->at(i_sci);
	const SCIDEIntoQConverter& cvt_first = par.qdc;
	
	if(cvt_first.IsDefaulted())
		ERROR("Tried to reach for the charge converter, but the initial file labelled '%s' has the defaulted converter?\n"
			"Full parameter is: %s\n", fileName.front().c_str(), mnd::streamable(par).c_str());

	/* Should be common for all the files.. If not, something is very bad. */
	FRSToFSingle const* const tof_p =
		   (i_sci == SCI21_I)                        ? tof_param->Get(SCI21_I, SCI22_I)
		: ((i_sci == SCI22_I) || (i_sci == SCI31_I)) ? tof_param->Get(SCI22_I, SCI31_I)
		:  (i_sci == SCI41_I)                        ? tof_param->Get(SCI31_I, SCI41_I)
		: nullptr;

	if(!tof_p)
		WARN("ToF parameter is nullptr, that's fine. Will not do β-correction in this case.");

	std::unique_ptr<phy::Nucleus> sec_nuc = beam_param->SecondaryIon();
	WARN("Secondary nucleus identified as: %s\n",
		sec_nuc->chem_to_string().c_str());

	const double brho = beam_param->BrhoAt(i_sci);
	double beta0 = phy::Beta(
		*sec_nuc, phy::Brho_t{brho}
	);
	const A2 beta_interval = {beta0 - beta_range, beta0 + beta_range};

	auto h1_sci_l = TH1P{Form("SCI%s QDC l [QDC units]", label[i_sci]),   0x009B2F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r = TH1P{Form("SCI%s QDC r [QDC units]", label[i_sci]),   0x008B1F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a = TH1P{Form("SCI%s E1 [QDC units]@Left-right combined, with #beta correction", label[i_sci]), 0x007B0F_c, bin[0], bin[1], bin[2]};

	auto h2_de0_vs_beta_all = TH2P(
		Form("((no_cut))E0 initial QDC [QDC units]:#beta [0..1]@SCI%s, no cut", label[i_sci]),
		beta_binning, beta_interval[0], beta_interval[1],
		bin[0], bin[1], bin[2]
	);
	auto h2_de1_vs_beta_all = TH2P(
		Form("((no_cut))E1 corrected QDC [QDC units]:#beta [0..1]@SCI%s, no cut", label[i_sci]),
		beta_binning, beta_interval[0], beta_interval[1],
		bin[0], bin[1], bin[2]
	);
	auto h2_de0_vs_beta = TH2P(
		Form("E0 initial QDC [QDC units]:#beta [0..1]@SCI%s", label[i_sci]),
		beta_binning, beta_interval[0], beta_interval[1],
		bin[0], bin[1], bin[2]
	);
	auto h2_de1_vs_beta = TH2P(
		Form("E1 corrected QDC [QDC units]:#beta [0..1]@SCI%s", label[i_sci]),
		beta_binning, beta_interval[0], beta_interval[1],
		bin[0], bin[1], bin[2]
	);
	auto h1_de0 = TH1P(
		Form("E0 initial QDC [QDC units]@SCI%s", label[i_sci]),
		mnd::col::rand_col(), bin[0], bin[1], bin[2]
	);
	auto h1_de1 = TH1P(
		Form("E1 corrected QDC [QDC units]@SCI%s", label[i_sci]),
		mnd::col::rand_col(), bin[0], bin[1], bin[2]
	);
	auto h1_foot_raw = TH1P{"((h1_raw))FOOT Q [charge]@No cut", 0xF27100_c, 140, 0, 7};
	auto h1_foot_cut = TH1P{"((h1_cut))FOOT Q [charge]@W/ cut", 0xB08100_c, 140, 0, 7};

	auto h1_sci_l_c = TH1P{Form("((h1_cut))SCI%s QDC l [QDC units]@With FOOT selection", label[i_sci]),   0x009B4F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r_c = TH1P{Form("((h1_cut))SCI%s QDC r [QDC units]@With FOOT selection", label[i_sci]),   0x007B2F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a_c = TH1P{Form("((h1_cut))SCI%s QDC avg [QDC units]@With FOOT selection", label[i_sci]), 0x005B0F_c, bin[0], bin[1], bin[2]};

	WARN("SCI%s: ToF parameter identified as: %s\n", label[i_sci],
		tof_p ? mnd::streamable(*tof_p).c_str() : "nullptr");
	
	ROOT::EnableThreadSafety();
	TApplication rootApp("app", 0, 0);

	show_console_cursor(false);
	DynamicProgress<ProgressBar> progress;

	mnd::parallel_process(fileName, nthreads, [=, &progress](size_t i, auto fname) mutable {
		
		/* Just to not have threads racing over mutable elements,
		 * we open the file here and fetch the parameter also on a per-file basis. */
		std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params_local;
		{
			std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.front().c_str(), "READ");
			get_obj(f, sci_params_local, "FRS_sci_parameters");
		}
		const SCIParam& p = sci_params_local->at(i_sci);
		const SCIDEIntoQConverter& cvt = p.qdc;
		
		/* Compare the two converters. They ought to be the same.. */
		if(cvt.IsDefaulted() || cvt != cvt_first)
			ERROR("Inside file: '%s', fetched a converter. \n"
				"      It is: %s\n"
				"Initial cvt: %s\n"
				"They are either not equal cvt, or it is invalid. Either way, not allowed!",
				fname.c_str(), mnd::streamable(cvt).c_str(), mnd::streamable(cvt_first).c_str());

		auto model = RNTupleModel::Create();

		std::shared_ptr<RNFOOTHit> foot {};
		if(q.is_some()) {
			foot = model->MakeField<RNFOOTHit>("FOOT");
		}
		auto frs = model->MakeField<RNFRSHit>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);
		const size_t nentries = ntuple->GetNEntries();

		auto _bar_ptr = std::make_unique<ProgressBar> (
			option::BarWidth{55},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("%zu/%zu: %zu (%s)", i+1, fileName.size(), nentries, fname.c_str())},
			option::ForegroundColor{ indicators::next_col() },
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];
		
		for(size_t entryId{0}; entryId < nentries; ++entryId ) {
			if( mnd::PrintProgress(bar, entryId, nentries, 1000) )
				progress.print_progress();

			ntuple->LoadEntry(entryId);
			const RNSciCal& sci = frs->cal.sci[i_sci];

			/* Calculate possible β */
			const double beta = (tof_p) ? (tof_p->Beta(frs->cal)) : NAN;

			const double e0_l  = p.E0_l(sci);
			const double e0_r  = p.E0_r(sci);
			const double e0    = p.E0(sci);
			const double e1    = p.E(sci,beta);

			h1_sci_l.Fill(e0_l);
			h1_sci_r.Fill(e0_l);
			h1_sci_a.Fill(e1);

			h2_de0_vs_beta_all.Fill(beta, e0);
			h2_de1_vs_beta_all.Fill(beta, e1);
		
			if(sci.hits.size() != 1)
				continue;

			if(q.is_some() && foot) {
				bool valid = false;
				
				for(const auto& track : foot->track) {
					h1_foot_raw->Fill(track.Q);
					if( mnd::IsInside(track.Q, qcut.unwrap()) ) {
						valid = true;
						h1_foot_cut->Fill(track.Q);
					}
				}
				if(!valid) continue;
			
				h1_sci_l_c->Fill(e0_l);
				h1_sci_r_c->Fill(e0_r);
				h1_sci_a_c->Fill(e1);
			}
			
			h2_de0_vs_beta->Fill(beta, e0);
			h2_de1_vs_beta->Fill(beta, e1);
			h1_de0->Fill(e0);
			h1_de1->Fill(e1);
		}
	});
	show_console_cursor(true);

	auto c = TCanvas{Form("sci%s",label[i_sci]), Form("SCI%s qdc values", label[i_sci]), 2180, 1420};
	c.Divide(4,3);
	c.cd( 1); h1_sci_l.Draw();
	c.cd( 2); h1_sci_r.Draw();
	c.cd( 3); h1_sci_a.Draw();
	c.cd( 4); h2_de0_vs_beta_all.Draw("COLZ"); gPad->SetLogz();
	c.cd( 8); h2_de1_vs_beta_all.Draw("COLZ"); gPad->SetLogz();
	c.cd( 9); h2_de0_vs_beta.Draw("COLZ");
	c.cd(10); h2_de1_vs_beta.Draw("COLZ");
	c.cd(11); h1_de0.Draw("HIST");
	c.cd(12); h1_de1.Draw("HIST");
	
	TCanvas* cf;
	if(q.is_some()) {
		cf = new TCanvas{"w_foot_cut", Form("SCI%s qdc values", label[i_sci]), 2000, 1300};
		cf->Divide(3,2);
		cf->cd(1); h1_sci_l_c.Draw();
		cf->cd(2); h1_sci_r_c.Draw();
		cf->cd(3); h1_sci_a_c.Draw();
		cf->cd(4); h1_foot_raw.Draw();
		cf->cd(5); h1_foot_cut.Draw();
	}

	WARN("End-of-main\n");
	rootApp.Run(); return 0;
}

