#include "util/CLI.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"

#include "TApplication.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

using MaybePed = mnd::Option<A2>;

int main(int argc, char* argv[]) {
    CLI::App app{"Calibrate the QDC into-charge measurement of SCI21/22/31. "
        "We don't do the velocity dependence correction, as it is anyway close to the minimum "
        "of the Bethe-Bloch curve. Maybe a TODO for later,..."};

    std::string fileName{};
	u32 i_sci = 0;
	A3 binning_sci = {4000, 0, 4000};
    u32 niter = 2;
    unsigned short line_size = 4;
    double sratio = 1.4;
    MaybePed ped = MaybePed::No;
    bool ped_from_file = false;
	auto save = canvas::Extension::nil;

    add_logged_option(app, "-f,--file", fileName, "Pass a file name.")
        ->check(CLI::ReadPermissions);
    add_logged_option(app, "-i,--sci", i_sci, 
        mnd::msg("Scintillator index; %u => SCI21, %u => SCI22, %u => SCI31",
            RNFRSCal::SCI21_I, RNFRSCal::SCI22_I, RNFRSCal::SCI31_I))
        ->check(CLI::Range(0, (int)RNFRSCal::SCI31_I));

	add_logged_option(app, "-s,--binning-sci", binning_sci, "Binning Y. If difference toggle given, then Y becomes the difference axis (delta axis).")
		->delimiter(',');
    add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder. Only with --ped option active.")
        ->check(CLI::PositiveNumber);
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak. Only with --ped option active.")
        ->check(CLI::PositiveNumber);
	add_logged_option(app, "-l,--line-size", line_size, "Fit curve line size.");
    add_logged_option(app, "-p,--ped", ped, "Do pedestal subtraction. Two numbers represent average pedestals for left and right channel.");
    add_logged_flag(app, "--ped-from-file", ped_from_file, "Take the pedestal values from the file itself. Will invalidate the --ped option's entered values.");
    add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileName.empty()) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}
    
    const auto& label = RNFRSCal::sci_label;
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
    {
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
        get_obj(f, sci_params, "FRS_sci_parameters");
    }
    [[maybe_unused]] const SCIParam& par = sci_params->at(i_sci);
    if( ped_from_file ) {
        /* In this case just take the pedestal values from the file. */
        u32 n_matched = par.SetConverter(fileName);
        if(n_matched == 0)
            ERROR("Requsted parameter from file \'%s\' itself, but no converter regex matched the file stem?\n",
                fileName.c_str());
        const SCIDEIntoQConverter* cvt = par.GetConverter();
        if(cvt->pedestal.left == 0 || cvt->pedestal.right == 0)
            ERROR("Requsted parameter from file \'%s\' itself, but parameter is left defaulted? Run without the flag first.\n",
                fileName.c_str());
        ped = MaybePed::Yes {
            cvt->pedestal.left,
            cvt->pedestal.right
        };
    }
    auto* h1_sci_ped_l = new TH1P (
        Form("((h1_sci))SCI%s-l pedestal [QDC units]", label[i_sci]),0xCB00CB_c, 2000, 0, 2000
    );
    auto* h1_sci_ped_r = new TH1P (
        Form("((h1_sci))SCI%s-r pedestal [QDC units]", label[i_sci]),0x13973F_c, 2000, 0, 2000
    );
    TH1P *h1_sci_l = nullptr, *h1_sci_r = nullptr, *h1_sci_e = nullptr;
    if( ped.is_some() ) {
        h1_sci_l = new TH1P (
            Form("((h1_sci))SCI%s-l value [QDC units]@Single hit cut", label[i_sci]), 0xABABAB_c, 4096, 0, 4096
        );
        h1_sci_r = new TH1P (
            Form("((h1_sci))SCI%s-r value [QDC units]@Single hit cut", label[i_sci]), 0xBABABA_c, 4096, 0, 4096
        );
        h1_sci_e = new TH1P (
            Form("((h1_sci))SCI%s average value [QDC units]@Single hit cut", label[i_sci]), 0xBABABA_c, 4096, 0, 4096
        );
    }
	TApplication rootApp("app", 0, 0);

    auto model = RNTupleModel::Create();
    auto frs = model->MakeField<RNFRSCal>("FRS");
    auto ntuple = (
        mnd::set_current_input_file(fileName),
        RNTupleReader::Open(std::move(model), "h103", fileName)
    );
    ProgressBar bar {
        option::BarWidth{50},
            option::Start{"["},
            option::Fill{"="},
            option::Lead{">"},
            option::Remainder{" "},
            option::End{"]"},
            option::PostfixText{mnd::msg("SCI%s QDC-cal", label[i_sci])},
            option::ForegroundColor{Color::yellow},
            option::ShowPercentage{true},
            option::ShowElapsedTime{true},
            option::ShowRemainingTime{true},
            option::FontStyles{std::vector{FontStyle::bold}}
    };
    const size_t nentries = ntuple->GetNEntries();

    for(auto entryId : *ntuple) {
        ntuple->LoadEntry(entryId);
        mnd::PrintProgress(bar, entryId, nentries, 500, mnd::dancer1, 0.30);
        
        const auto& sci = frs->sci[i_sci];
        
        if(sci.hits.empty()) {
            h1_sci_ped_l->Fill(sci.El);
            h1_sci_ped_r->Fill(sci.Er);
            continue;
        }
        
        if(ped.is_some()) {
            /* Demand single hit entries. */
            if(sci.hits.size() != 1) continue;

            const f64 de_l = std::max(sci.El - ped.unwrap()[0], 0.0);
            const f64 de_r = std::max(sci.Er - ped.unwrap()[1], 0.0);
            h1_sci_l->Fill(de_l);
            h1_sci_r->Fill(de_r);
            h1_sci_e->Fill( sqrt(de_l * de_r) );
        }
    }

    /* In either case, just try to match the pedestal first. */
    TCanvas *cp = new TCanvas("Pedestal", "Pedestal", 1800, 1200);
    cp->Divide(2,1);
    cp->cd(1);
    auto [fit_l, _err_l] = h1_sci_ped_l->DrawAndFit(sratio, kRed, line_size, niter);
    cp->cd(2);
    auto [fit_r, _err_r] = h1_sci_ped_r->DrawAndFit(sratio, kRed, line_size, niter);

    SCIQDCPedestal result;
    result.left  = fit_l[1];
    result.right = fit_r[1];
    WARN("Pedestal of SCI%s found:\n", label[i_sci]);
    std::cerr << "\"pedestal\": " << nlohmann::json(result).dump(4) << std::endl;

    if(ped.is_some()) {
        TCanvas *c = new TCanvas("qdc_sub", "qdc_sub", 1800, 1200);
        c->Divide(2,2);
        c->cd(1);
        h1_sci_l->Draw();
        c->cd(2);
        h1_sci_r->Draw();
        c->cd(3);
        h1_sci_e->Draw();
    }

    WARN("End-of-main");
    rootApp.Run(); return 0;
}
