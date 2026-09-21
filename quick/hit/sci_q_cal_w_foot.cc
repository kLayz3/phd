#include "TFRSCalCont.h"
#include "util/CLI.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"

#include "TApplication.h"
#include "TFRSHitCont.h"
#include "TFOOTHitCont.h"
#include <algorithm>

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

/* I'm tilted and writing this script like an ape. Entschuldigung for my caveman behaviour. */
int main(int argc, char* argv[]) {
    CLI::App app{"Calibrate the QDC into-charge measurement of SCI21/22/31. "
        "We use FOOT's tracking here to imply the charge Q."};

	std::vector<std::string> fileName;
	u32 nthreads = 1;
	u32 q = 6;
	double d = 0.4;
	A2 ped21 = {NAN,NAN};
	A2 ped22 = {NAN,NAN};
	A2 ped31 = {NAN,NAN};
	A3 bin = {400,100,4000};

	add_logged_option(app, "-f,--file", fileName, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-q,--charge", q, "Select the charge: 1...=6 on which to cut on FOOT's")
		->check(CLI::Range(1,6));
	add_logged_option(app, "-d,--range", d, "Select the charge interval width, [q-d, q+d], on which to cut on FOOT's")
		->check(CLI::Range(0.01,0.50));
	add_logged_option(app, "-b,--binning", bin, "Binning for SCI QDC axis.")
		->delimiter(',');
    add_logged_option<DisplayDefault::No>(app, "--ped21", ped21, "Two numbers represent average pedestals for left and right channel, respectively of SCI21.")
		->delimiter(',')
		->check(CLI::Range(1.0, 200.0))
		->required();
    add_logged_option<DisplayDefault::No>(app, "--ped22", ped22, "Two numbers represent average pedestals for left and right channel, respectively of SCI22.")
		->delimiter(',')
		->check(CLI::Range(1.0, 200.0))
		->required();
    add_logged_option<DisplayDefault::No>(app, "--ped31", ped31, "Two numbers represent average pedestals for left and right channel, respectively of SCI31.")
		->delimiter(',')
		->check(CLI::Range(1.0, 200.0))
		->required();

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;
	if(fileName.size() == 0)
		ERROR("To continue, must supply at least one file name!\n");
	
	const A2 qcut = std::array{q-d, q+d};

	ROOT::EnableThreadSafety();
	TApplication rootApp("app", 0, 0);

	auto h1_sci_l21 = TH1P{"SCI21 QDC r [QDC units]",   0xCC00CC_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r21 = TH1P{"SCI21 QDC l [QDC units]",   0xBC00BB_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a21 = TH1P{"SCI21 QDC avg [QDC units]", 0xAC00AB_c, bin[0], bin[1], bin[2]};
	auto h1_sci_l22 = TH1P{"SCI22 QDC l [QDC units]",   0x0070DD_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r22 = TH1P{"SCI22 QDC l [QDC units]",   0x0060CD_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a22 = TH1P{"SCI22 QDC avg [QDC units]", 0x0050BD_c, bin[0], bin[1], bin[2]};
	auto h1_sci_l31 = TH1P{"SCI31 QDC l [QDC units]",   0x009B2F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r31 = TH1P{"SCI31 QDC r [QDC units]",   0x008B1F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a31 = TH1P{"SCI31 QDC avg [QDC units]", 0x007B0F_c, bin[0], bin[1], bin[2]};

	auto h1_sci_l_c21 = TH1P{"((h1_cut))SCI21 QDC r [QDC units]@ w cut",   0xCC00CC_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r_c21 = TH1P{"((h1_cut))SCI21 QDC l [QDC units]@ w cut",   0xAC00AB_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a_c21 = TH1P{"((h1_cut))SCI21 QDC avg [QDC units]@ w cut", 0x7C008B_c, bin[0], bin[1], bin[2]};
	auto h1_sci_l_c22 = TH1P{"((h1_cut))SCI22 QDC l [QDC units]@ w cut",   0x0070DD_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r_c22 = TH1P{"((h1_cut))SCI22 QDC l [QDC units]@ w cut",   0x0050BD_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a_c22 = TH1P{"((h1_cut))SCI22 QDC avg [QDC units]@ w cut", 0x00308D_c, bin[0], bin[1], bin[2]};
	auto h1_sci_l_c31 = TH1P{"((h1_cut))SCI31 QDC l [QDC units]@ w cut",   0x009B4F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_r_c31 = TH1P{"((h1_cut))SCI31 QDC r [QDC units]@ w cut",   0x007B2F_c, bin[0], bin[1], bin[2]};
	auto h1_sci_a_c31 = TH1P{"((h1_cut))SCI31 QDC avg [QDC units]@ w cut", 0x005B0F_c, bin[0], bin[1], bin[2]};

	auto h1_foot_raw = TH1P{"((h1_raw))FOOT Q [charge]@No cut", 0xF27100_c, 140, 0, 7};
	auto h1_foot_cut = TH1P{"((h1_cut))FOOT Q [charge]@W/ cut", 0xB08100_c, 140, 0, 7};

	show_console_cursor(false);
	DynamicProgress<ProgressBar> progress;

	mnd::parallel_process(fileName, nthreads, [=, &progress](size_t i, auto fname) mutable {
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTHit>("FOOT");
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
			const auto& sci21 = frs->cal.sci[0];
			const auto& sci22 = frs->cal.sci[1];
			const auto& sci31 = frs->cal.sci[2];

#define GET_N_FILL_LETSGO(NUM) \
            const f64 de_l_##NUM = std::max(sci##NUM.El - ped##NUM[0], 0.0); \
            const f64 de_r_##NUM = std::max(sci##NUM.Er - ped##NUM[1], 0.0); \
			h1_sci_l##NUM->Fill(de_l_##NUM); \
			h1_sci_r##NUM->Fill(de_r_##NUM); \
			h1_sci_a##NUM->Fill( sqrt(de_l_##NUM * de_r_##NUM) );
			
#define FILL_CUT_LETSGO(NUM) \
			if(sci##NUM.hits.size() == 1) { \
				h1_sci_l_c##NUM->Fill(de_l_##NUM); \
				h1_sci_r_c##NUM->Fill(de_r_##NUM); \
				h1_sci_a_c##NUM->Fill( sqrt(de_l_##NUM * de_r_##NUM) ); \
			}

			GET_N_FILL_LETSGO(21)
			GET_N_FILL_LETSGO(22)
			GET_N_FILL_LETSGO(31)
			
			bool valid {false};
			for(const auto& track : foot->track) {
				h1_foot_raw->Fill(track.Q);
				if( mnd::IsInside(track.Q, qcut) ) {
					valid = true;
					h1_foot_cut->Fill(track.Q);
				}
			}
			if(!valid) continue;
			
			FILL_CUT_LETSGO(21)
			FILL_CUT_LETSGO(22)
			FILL_CUT_LETSGO(31)
		}
	});
	show_console_cursor(true);

#define DECL_CANVAS_LETSGO(NUM) \
	TCanvas *c##NUM = new TCanvas("SCI" #NUM , "SCI" #NUM " qdc values before and after FOOT cut.", 2150, 1400); \
	c##NUM->Divide(3,2); \
	c##NUM->cd(1); h1_sci_l##NUM.Draw(); \
	c##NUM->cd(2); h1_sci_r##NUM.Draw(); \
	c##NUM->cd(3); h1_sci_a##NUM.Draw(); \
	c##NUM->cd(4); h1_sci_l_c##NUM.Draw(); \
	c##NUM->cd(5); h1_sci_r_c##NUM.Draw(); \
	c##NUM->cd(6); h1_sci_a_c##NUM.Draw();

	DECL_CANVAS_LETSGO(21)
	DECL_CANVAS_LETSGO(22)
	DECL_CANVAS_LETSGO(31)

	TCanvas *cf = new TCanvas("FOOT", "FOOT", 2000, 1300);
	cf->Divide(2,1);
	cf->cd(1); h1_foot_raw->Draw();
	cf->cd(2); h1_foot_cut->Draw();

	WARN("End-of-main\n");
	rootApp.Run(); return 0;
}

