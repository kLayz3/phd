#include "util/CLI.h"

#include "util/MacroHelpers.h"
#include "util/PrettyHisto.hxx"

#include "TApplication.h"
#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

struct Dead { A2 val; };
std::istream& operator>>(std::istream& , Dead& );
std::ostream& operator<<(std::ostream& , const Dead& );

int main(int argc, char* argv[]) {
	CLI::App app{"Calibrate referent FOOT detectors' offset, based on the primary beam cal run.\n\
	              This alignment algorithm is reserved only for non-rotated detectors.\n\
                  For referent detectors, use the `fit_angle` routine."};
	
	std::vector<std::string> fileName;
	uint32_t ifoot = -1;
	std::array<double,3> binning_x = {200,-30,30};
	std::array<double,2> foot_q_cut = {5.4, 6.6};
	std::vector<Dead> dead_regions = {};
	std::array<double,2> sci21_cut = {NAN, NAN};
	std::array<double,2> sci22_cut = {NAN, NAN};
	std::array<double,2> sci31_cut = {NAN, NAN};
	auto save = canvas::Extension::nil;
	bool do_verify = false;
	size_t max_events = -1; 

	add_logged_option(app, "-f,--file", fileName, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-i,--foot-id", ifoot, 
		"Select which FOOT detector.")
		->check([](const std::string& match) -> std::string { 
			u32 id = -1;
			mnd::parse(match, id);
			return (id > 7 or (id == 4) or (id == 5))? "Must be {0,1,2,3, 6,7}": "";
		})
		->mandatory();
	add_logged_option<DisplayDefault::No>(app, "-m,--max-events", max_events, 
		"Max events (at most) taken from each ROOT file. Default: all entries.");

	add_logged_option<DisplayDefault::No>(app, "-d,--dead", dead_regions,
		"Dead regions [in mm] that won't go into the fit, separated by \';\'")
		->delimiter(';')
		->type_name("LO,HI[;LO,HI...]");
	add_logged_option(app, "-b,--bins-x",binning_x, "Binning X")
		->delimiter(',');
	add_logged_option(app, "-q,--foot-cut",foot_q_cut, "FOOT Q cut (charge)")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(','); 
	add_logged_flag(app, "--verify", do_verify, "Already take inputted parameters and calculate the true position.");
	add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(fileName.size() == 0) {
		WARN("To continue, must supply at least one file name!\n"); return 0;
	}

	TApplication rootApp("app", 0, 0);
	FOOTParam *foot_param; 
	FOOTBoxParam *box;
	{
		const auto& fname = fileName.front();
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fname.c_str(), "READ");
		get_obj(f, foot_param, Form("FOOT%d_setup", ifoot));
		get_obj(f, box, "FOOT0_box");
	}
	Orientation o = foot_param->GetOrientation();
	if(o == Orientation::UNKNOWN)
		ERROR("FOOT orientation not specified. I won't allow it.\n");

	auto* foot_q_vs_d = new TH2P(Form("Cluster Charge:Delta [-0.5, 0.5]@FOOT%d", ifoot), 
		80, -0.5, 0.5, 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_q_vs_x = new TH2P(Form("Cluster Charge:FOOT measurement [mm]@FOOT%d", ifoot), 
		binning_x[0], binning_x[1], binning_x[2], 100, foot_q_cut[0], foot_q_cut[1]);
	auto* foot_pos = new TH1P(Form("((h1))FOOT measurement [mm]@FOOT%d", ifoot), kCyan - 6, binning_x[0], binning_x[1], binning_x[2]);
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", ORGB{0xCB00CB}, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", ORGB{0x0070DD}, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", ORGB{0x009B2F}, 500, 300, 4000);
	auto* h1_sci21_cut = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", ORGB{0x890389}, 500, 300, 4000);
	auto* h1_sci22_cut = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", ORGB{0x6180FD}, 500, 300, 4000);
	auto* h1_sci31_cut = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", ORGB{0x7DE69D}, 500, 300, 4000);

	for(size_t i{0}; i < fileName.size(); ++i) {
		const auto& fname = fileName[i];
		auto model = RNTupleModel::Create();
		auto frs  = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
		auto foot = model->MakeField<RNFOOTCal>(Form("FOOT%d", ifoot));
		auto ntuple = RNTupleReader::Open(std::move(model), "h103", fname);

		const size_t nentries = ( (max_events < ntuple->GetNEntries()) ? max_events : ntuple->GetNEntries() );
	
		ProgressBar bar {
			option::BarWidth{50},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{"~"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("Analysis (per event: %s)", fname.c_str())},
			option::ForegroundColor{Color::green},
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		};
		WARN("Proceeding with file [%zu/%zu]: \'%s\'. Entries: [%'zu]\n", i+1, fileName.size(), fname.c_str(), nentries);
		
		for(size_t entryId{0}; entryId < nentries; ++entryId ) {
			mnd::PrintProgress(bar, entryId, nentries, 1000);

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

			if(sci21.hits.size() >= 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31_cut->Fill(sci31.E);

			/* In an event, only a single valid FOOT cluster must be found. */
			bool is_valid = false;
			double xFOOT_ = NAN;
			double qFOOT_ = NAN;
			double dFOOT_ = NAN;
			for(const auto& hit : foot->fCl) {
				if(hit.fCM == 1) continue; 	
				double q = foot_param->Q( hit );
				if(mnd::IsValid(foot_q_cut) and !mnd::IsInside(q, foot_q_cut)) continue;

				double d = hit.Delta();

				double hit_position;
				if(do_verify == true)
					hit_position = foot_param->X0(hit);
				else
					hit_position = foot_param->BarePosition(hit);

				if( std::isfinite(xFOOT_) ) {
					/* Already found valid point in the event. */
					is_valid = false; break;
				} else {
					/* Export it outside. */
					xFOOT_ = hit_position;
					qFOOT_ = q;
					dFOOT_ = d;
				}
				is_valid = true;
			}
			if(!is_valid) continue;

			foot_pos->Fill(xFOOT_);
			foot_q_vs_d->Fill(dFOOT_, qFOOT_);
			foot_q_vs_x->Fill(xFOOT_, qFOOT_);
		}
		bar.mark_as_completed();
	}
	show_console_cursor(true);

	/* Some strips/regions are dead and not to confuse the gauss-chan 🥺 👉👈,
	 * exclude these regions from the fit. Passed in as a region sequence. */
	constexpr double fit_area_factor = 0.5;
	auto mygaus = [dead_regions](Double_t* x, Double_t* par) -> Double_t {
		double xx = x[0];
		for(const auto& region : dead_regions) {
			if(mnd::IsInside(xx, region.val)) {
				TF1::RejectPoint();
				return 0.0;
			}
		}
		return par[0] * std::exp(
			-0.5 * std::pow((xx-par[1])/par[2], 2)
		);
	};

	TF1* f = new TF1("f", mygaus, fit_area_factor*binning_x[1], fit_area_factor*binning_x[2], 3);
	f->SetParNames("A", "mu", "sigma");
    f->SetParameters((*foot_pos)->GetMaximum(), (*foot_pos)->GetMean(), (*foot_pos)->GetStdDev());
	f->SetBit(TObject::kCanDelete, false);
	f->SetBit(TF1::kNotGlobal);

	(*foot_pos)->Fit(f, "RMN");
	WARN("Result to be put in setup file: " BOLD "\"delta_p\": %.5f\n" KNRM, -f->GetParameter("mu"));
	TCanvas* c = new TCanvas("FOOT", Form("FOOT%d calibration position plots", ifoot), 2150,1400);
	c->Divide(2,2);
	c->cd(1); foot_q_vs_d->Draw("COLZ"); gPad->SetLogz();
	c->cd(2); foot_q_vs_x->Draw("COLZ"); gPad->SetLogz();
	c->cd(3); foot_pos->Draw("COLZ");
	
	f->Draw("same");

	TCanvas* cs = new TCanvas("SCIs&TPCs", "SCI21,22,31 and TPC ref", 2200, 1200);
	cs->Divide(3,2);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();

	canvas::save_all<canvas::Exe>(save, { Form("FOOT%d", ifoot) });

	WARN("End-of-main");
	rootApp.Run(); return 0;
}

std::istream& operator>>(std::istream& in, Dead& out) {
	return ::mnd::template operator>> <','>(in, out.val);
}
std::ostream& operator<<(std::ostream& os, const Dead& out) {
	return os << out.val;
}
