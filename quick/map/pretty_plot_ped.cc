#include "util/CLI.h"
#include "monad/monad.hxx"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"
#include "util/TextUtils.h"

#include "TApplication.h"
#include "TGraph.h"

#include "TFOOTMapCont.h"
#include "TFRSMapCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

constexpr auto N_ASIC   = RNFOOTMap::N_ASIC;
constexpr auto N_STRIPS = RNFOOTMap::N_STRIPS;
constexpr auto N_STRIPS_PER_ASIC = RNFOOTMap::N_STRIPS_PER_ASIC;

int main(int argc, char* argv[]) {
	ROOT::EnableThreadSafety();

	CLI::App app{"Pretty plot raw FOOT spectrum, with and without pedestal subtraction."};

	std::vector<std::string> fileNames {};
	i32 N = 0;

	A3 foot_binning = {3000, -500, 2500};
	A3 asic_offset_binning = {5000,-2500,2500};

	A2 sci21_cut = {NAN, NAN};
	A2 sci22_cut = {NAN, NAN};
	A2 sci31_cut = {NAN, NAN};
	A3 sci_binning = {500, 300, 4000};
	
	std::vector<canvas::Extension> save{};
	u32 nthreads = 1;
	
	add_logged_option(app, "-f,--file", fileNames, "Pass a file name.")
		->check(CLI::ReadPermissions)
        ->delimiter(',');
	add_logged_option(app, "-i,--de10-id", N,
		"Select which FOOT detector, based on the DE10 label.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "--sci-binning",sci_binning, "Scintillator QDC binning")
		->delimiter(',');
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
	)->check(CLI::PositiveNumber);
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileNames.size() == 0) {
		WARN("To continue, must supply a valid file name!\n"); return 0;
	}

	using SigmaType = mnd::StatisticalPoint<double, u64>;
	std::array<
		std::vector<SigmaType>, N_STRIPS
	> ped_sigma {};
	std::array<double, N_STRIPS>* tmp_csig;

	TH2I* raw = nullptr; TH2D* corr = nullptr;
	for(const auto& fname : fileNames) {
		auto f = std::make_unique<TFile>(fname.c_str(), "READ");
		if(!raw) {
			get_obj(f, raw,  Form("FOOT%d_h2_raw",  N));
		} else {
			TH2I* tmp;
			get_obj(f, tmp,  Form("FOOT%d_h2_raw",  N));
			raw->Add(tmp);
		}
		if(!corr) {
			get_obj(f, corr,  Form("FOOT%d_h2_corr",  N));
		} else {
			TH2D* tmp;
			get_obj(f, tmp,  Form("FOOT%d_h2_corr",  N));
			corr->Add(tmp);
		}

		get_obj(f, tmp_csig, Form("FOOT%d_ped_sigma", N));

		const auto ntuple = RNTupleReader::Open(RNTupleModel::Create(), "h104", fname);
		const u64 nentries = ntuple->GetNEntries();

		for(int i=0; i<N_STRIPS; ++i) {
			auto& values = ped_sigma[i];
			values.emplace_back( SigmaType{
				.value  = tmp_csig->operator[](i),
				.weight = nentries
			});
		}
	}
	TApplication rootApp("app", 0, 0);
	
	auto h1_sci21 = TH1P{"SCI21 QDC mean [QDC units]", 0x6180FD_c, sci_binning};
	auto h1_sci22 = TH1P{"SCI22 QDC mean [QDC units]", 0x0070DD_c, sci_binning};
	auto h1_sci31 = TH1P{"SCI31 QDC mean [QDC units]", 0x009B2F_c, sci_binning};
	auto h1_sci21_cut = TH1P{"((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, sci_binning};
	auto h1_sci22_cut = TH1P{"((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0xCB00CB_c, sci_binning};
	auto h1_sci31_cut = TH1P{"((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, sci_binning};

	auto corr_cut = TH2P{"((h2_cut))FOOT E [ADC units]:Strip mumber@With sci cut",
		N_STRIPS, 0, N_STRIPS, foot_binning};
	auto h2_ped_offset = TH2P{"((h2_offset))Common ASIC offset [ADC units]:ASIC number [0...9]@With cut", 
		N_ASIC, 0, N_ASIC, asic_offset_binning};


	DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
        auto model = RNTupleModel::Create();
        auto foot = model->MakeField<RNFOOTMap>(Form("FOOT%d", N));
		auto frs = model->MakeField<RNFRSMap>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);

        const size_t nentries = ntuple->GetNEntries();

		auto _bar_ptr = std::make_unique<ProgressBar>(
            option::BarWidth{50},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PostfixText{mnd::msg(
					"%zu/%zu: %'zu (%s)", i+1, fileNames.size(), nentries, fname.c_str()
				)},
                option::ForegroundColor{ indicators::next_col() },
                option::ShowPercentage{true},
                option::ShowElapsedTime{true},
                option::ShowRemainingTime{true},
                option::FontStyles{std::vector{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];

        for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
            if(mnd::PrintProgress(bar, entryId, nentries, 1000))
				progress.print_progress();
			const auto& sci21 = frs->sci[0];
			const auto& sci22 = frs->sci[1];
			const auto& sci31 = frs->sci[2];
			
			const double e21 = sqrt((double)sci21.qdc[0] * sci21.qdc[1]);
			const double e22 = sqrt((double)sci22.qdc[0] * sci22.qdc[1]);
			const double e31 = sqrt((double)sci31.qdc[0] * sci31.qdc[1]);

			if(sci21.tdc.size() >= 1) h1_sci21->Fill(e21);
            if(sci22.tdc.size() >= 1) h1_sci22->Fill(e22);
            if(sci31.tdc.size() >= 1) h1_sci31->Fill(e31);

            if(mnd::IsValid(sci21_cut) and (sci21.tdc.size() != 1 or !mnd::IsInside(e21, sci21_cut))) continue;
            if(mnd::IsValid(sci22_cut) and (sci22.tdc.size() != 1 or !mnd::IsInside(e22, sci22_cut))) continue;
            if(mnd::IsValid(sci31_cut) and (sci31.tdc.size() != 1 or !mnd::IsInside(e31, sci31_cut))) continue;

			if(mnd::IsValid(sci21_cut) and (sci21.tdc.size() != 1 or !mnd::IsInside(e21, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.tdc.size() != 1 or !mnd::IsInside(e22, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.tdc.size() != 1 or !mnd::IsInside(e31, sci31_cut))) continue;

			if(sci21.tdc.size() == 1) h1_sci21_cut->Fill(e21);
			if(sci22.tdc.size() == 1) h1_sci22_cut->Fill(e22);
			if(sci31.tdc.size() == 1) h1_sci31_cut->Fill(e31);

			for(int i=0; i<N_STRIPS; ++i)
				corr_cut->Fill(i, foot->FOOTE[i]);
			for(int a=0; a<N_ASIC; ++a)
				h2_ped_offset->Fill(a, foot->common_offset[a]);
		}
		if(!bar.is_completed())
			bar.mark_as_completed();
    });

	std::vector<TLine*> vlines;
	for(int i = 1; i < 10; ++i) {
		TLine* line = new TLine(i * N_STRIPS_PER_ASIC, -500,
		                        i * N_STRIPS_PER_ASIC, 4096);
		line->SetLineColor(kRed);
		line->SetLineStyle(2);
		line->SetLineWidth(3);
		vlines.push_back( line );
	}

	TCanvas* cped;
	cped = new TCanvas("cped", "Pedestal", 1000, 2000);
	cped->Divide(3,1);

	cped->cd(1);
	gPad->SetLogz();
	raw->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(raw->GetYaxis()->GetXmin());
		l->SetY2(raw->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}

	cped->cd(2);
	gPad->SetLogz();
	corr->Draw("COLZ");
	for(auto* l0 : vlines) {
		TLine* l = dynamic_cast<TLine*>(l0->Clone());
		l->SetY1(corr->GetYaxis()->GetXmin());
		l->SetY2(corr->GetYaxis()->GetXmax());
		l->Draw("SAME");
	}
	cped->cd(3); gPad->SetLogz();
	corr_cut->Draw("COLZ");

	const std::array<double, N_STRIPS> final_sigma
		= mnd::map(ped_sigma, [](const auto& values) {
			return mnd::mean_var(values).mean;
		});

	TGraph* s0 = new TGraph(N_STRIPS, final_sigma.data());
	s0->SetMarkerStyle(20);
	s0->SetMarkerSize(1.1);
	s0->SetMarkerColor(kRed - 1);
	TCanvas* csig = new TCanvas("csig", "Sigma", 1618, 1000);
	csig->Divide(2,1);
	csig->cd(1);
	s0->SetTitle(Form("Pedestal Width: FOOT%d;Strip number;ADC value", N));
	s0->Draw("AP");
	csig->cd(2);
	h2_ped_offset->Draw("COLZ");
	
	TCanvas* cs = new TCanvas("cs", "SCI21,22,31", 2000, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	canvas::save_all<canvas::Exe>(save, {
		Form("FOOT%d", N),
		mnd::rootfile_names_concatenated(fileNames)
	});
	WARN("End-of-main\n");
	rootApp.Run();

    return 0;
}
