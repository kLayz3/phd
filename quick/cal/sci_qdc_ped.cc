#include "TFRSCalCont.h"
#include "nlohmann/json.hpp"
#include "util/CLI.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"

#include "TApplication.h"

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
    CLI::App app{"Find the QDC pedestals of SCI21/22/31/41. "};

	std::vector<std::string> fileNames{};
	u32 i_sci = 0;
    u32 niter = 2;
    unsigned short line_size = 4;
    double sratio = 1.4;
	auto save = canvas::Extension::nil;

	add_logged_option(app, "-f,--file", fileNames, "Pass a file name(s).")
        ->check(CLI::ReadPermissions)
		->delimiter(',');
    add_logged_option(app, "-i,--sci", i_sci,
        mnd::msg("Scintillator index; %u => SCI21, %u => SCI22, %u => SCI31, %u => SCI41",
            SCI21_I, SCI22_I, SCI31_I, SCI41_I))
        ->check(CLI::Range(
			(int)SCI21_I, (int)SCI41_I
		));

    add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder. Only with --ped option active.")
        ->check(CLI::PositiveNumber);
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak. Only with --ped option active.")
        ->check(CLI::PositiveNumber);
	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
    add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.empty()) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

    const auto& label = RNFRSCal::sci_label;
	
	auto* h1_sci_ped_l = new TH1P (
        Form("((h1_sci))SCI%s-l pedestal [QDC units]", label[i_sci]),0xCB00CB_c, 2000, 0, 2000
    );
    auto* h1_sci_ped_r = new TH1P (
        Form("((h1_sci))SCI%s-r pedestal [QDC units]", label[i_sci]),0x13973F_c, 2000, 0, 2000
    );


	auto* h1_sci_l = new TH1P (
		Form("((h1_sci))SCI%s-l value [QDC units]@Single hit cut", label[i_sci]), 0xABABAB_c,
			4096, 0, 4096);
	auto* h1_sci_r = new TH1P (
		Form("((h1_sci))SCI%s-r value [QDC units]@Single hit cut", label[i_sci]), 0xBABABA_c,
			4096, 0, 4096);
	auto* h1_sci_e = new TH1P (
		Form("((h1_sci))SCI%s average value [QDC units]@Single hit cut", label[i_sci]), 0xBABABA_c,
			4096, 0, 4096);

	TApplication rootApp("app", 0, 0);

	for(const auto& file : fileNames) {
		auto model = RNTupleModel::Create();
		auto frs = model->MakeField<RNFRSCal>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h103", file);
		ProgressBar bar {
			option::BarWidth{50},
				option::Start{"["},
				option::Fill{"="},
				option::Lead{">"},
				option::Remainder{" "},
				option::End{"]"},
				option::PostfixText{mnd::msg("SCI%s QDC-cal [%s]", label[i_sci], file.c_str())},
				option::ForegroundColor{Color::yellow},
				option::ShowPercentage{true},
				option::ShowElapsedTime{true},
				option::ShowRemainingTime{true},
				option::FontStyles{std::vector{FontStyle::bold}}
		};
		const size_t nentries = ntuple->GetNEntries();
	
		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			mnd::PrintProgress(bar, entryId, nentries, 500);
			
			const auto& sci = frs->sci[i_sci];

			if(sci.hits.empty()) {
				h1_sci_ped_l->Fill(sci.El);
				h1_sci_ped_r->Fill(sci.Er);
				continue;
			} else {
				h1_sci_l->Fill(sci.El);
				h1_sci_r->Fill(sci.Er);
				h1_sci_e->Fill(std::sqrt(sci.El * sci.Er) );
			}
		}
	}
    TCanvas *cp = new TCanvas("pedestal", "Pedestal", 1800, 1200);
    cp->Divide(2,1);
    cp->cd(1);
    auto [fit_l, _err_l] = h1_sci_ped_l->DrawAndFit(sratio, kRed, line_size, niter);
    cp->cd(2);
    auto [fit_r, _err_r] = h1_sci_ped_r->DrawAndFit(sratio, kRed, line_size, niter);

	SCIDEIntoQConverter qcvt {};
	SCIQDCPedestal& result = qcvt.pedestal;
    result.left  = fit_l[1];
    result.right = fit_r[1];
	std::string key = mnd::longest_common_prefix(fileNames);
	nlohmann::json j{};
	j[key] = result;
    WARN("Pedestal of SCI%s found:\n%s%s%s\n", label[i_sci],
		 BOLD, j.dump(4).c_str(), KNRM);

    SCIQDCPedestal peak;
	TH1D* hl = &h1_sci_ped_l->h;
	TH1D* hr = &h1_sci_ped_r->h;
	result.left  = hl->GetXaxis()->GetBinCenter( hl->GetMaximumBin() );
	result.right = hr->GetXaxis()->GetBinCenter( hr->GetMaximumBin() );
    WARN("Pedestal PEAK of SCI%s found:\n%s\n", label[i_sci], j.dump(4).c_str());

    TCanvas *c = new TCanvas("qdc-raw", "Raw QDC (for hits)", 1800, 1200);
	c->Divide(3,1);
	c->cd(1); h1_sci_l->Draw();
	c->cd(2); h1_sci_r->Draw();
	c->cd(3); h1_sci_e->Draw();

	canvas::save_all<canvas::Exe>( save, { mnd::msg("SCI%d", i_sci) });

    WARN("End-of-main\n");
    rootApp.Run(); return 0;
}
