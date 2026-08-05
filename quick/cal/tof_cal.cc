/* Calibrate the time-of-flight between: 
 * [0] : SCI21 -> SCI22 , to indicate "possible" velocity at S2 
 * [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity at S3. 
 *
 * No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run. 
 * All of these measurements must be done with 12C files. There are 4 of relevance. */

#include "util/CLI.h"

#include "util/GaussFitMax.hxx"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.hxx"
#include "IonOptics.hxx"

#include "TApplication.h"
#include "TFRSCalCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

struct Brho {
    static constexpr char SEP = ';';
    double s2_incoming() const noexcept {
        return data_[0];
    }
    double s2_outgoing() const noexcept {
        return data_[1];
    };
    std::array<double, 2> data_;
};

struct FileBrho {
    static constexpr char SEP = ':';
    std::string name;
    Brho brho;
};
std::istream& operator>>(std::istream& , Brho& );
std::ostream& operator<<(std::ostream& , const Brho& );
std::istream& operator>>(std::istream& , FileBrho& );
std::ostream& operator<<(std::ostream& , const FileBrho& );

 int main(int argc, char* argv[]) {
	CLI::App app{"Calibrate the time-of-flight between:\n\
                  [0] : SCI21 -> SCI22 , to indicate \"possible\" velocity at S2\n\
                  [1] : SCI22 -> SCI31 , optional, due to low transmission, to measure velocity at S3.\n\
                  [2] : SCI21 -> SCI31 , as a intermediate calculation for [0].\n\
                  No need to look at 21-41 or 21-31 ToF, as this ToF is anyway invalid for main 9C run.\n\
                  All of these measurements must be done with 12C files. There are 3 or 4 files of relevance."};
    
    
	std::vector<FileBrho> f;
	std::array<double,3> dt_cut21_31 = {1000, -100, 100};
	std::array<double,3> dt_cut22_31 = {1000, -100, 100};
    double sratio = GAUSS_FIT_SIDE_RATIO_DEFAULT;
    double niter = 2;
	auto save = canvas::Extension::nil;
    add_logged_option(app, "-f,--file", f, "Pass one or more file names and corresponding 2 brho's.")
		->delimiter(',')
        ->type_name("[NAME:BRHO1;BRHO2 , ...]");

	add_logged_option<DisplayDefault::No>(app, "--dt-cut-21-31", dt_cut21_31, 
        "Delta T (SCI31-SCI21) cut, in TDC units [25ps]")
		->delimiter(','); 
	add_logged_option<DisplayDefault::No>(app, "--dt-cut-22-31", dt_cut22_31, 
        "Delta T (SCI31-SCI22) cut, in TDC units [25ps]")
		->delimiter(','); 
    add_logged_option(app, "--sratio", sratio, "Width ratio of raw histogram, how much to fit around the peak.")
		->check(CLI::PositiveNumber); 
	add_logged_option(app, "--niter", niter, "Gaussian TH1D fit, number of iterations for the peak finder.")
		->check(CLI::PositiveNumber);
	add_enum_option(app, "-o,--save", save, "Save the resulting histogram as an extension.");

    bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

	if(f.size() < 3) {
		WARN("To continue, must supply at least 3 file names!\n"); return 0;
	}
    
	TApplication rootApp("app", 0, 0);

    const u32 SCI_21_I = TFRSCalCont::sci_moniker.at("21");
    const u32 SCI_22_I = TFRSCalCont::sci_moniker.at("22");
    const u32 SCI_31_I = TFRSCalCont::sci_moniker.at("31");
    std::vector <
        std::pair<TH1P*, TH1P*>
    > hist; 

	std::vector<double> x, y; // fitting containers
    u32 i = 0;
	for(const auto& [fileName, brho] : f) {
        ++i;
        auto model = RNTupleModel::Create();
        auto frs = model->MakeField<RNFRSCal>("FRS");
        auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
        const double beta_inc = phy::Beta(6, 12, brho.s2_incoming());
        const double beta_out = phy::Beta(6, 12, brho.s2_outgoing());

		auto* h1_dt21_31 = new TH1P(Form("((h1_dt%d_21_31))Delta t [25 ps]@SCI31 - SCI21@TOF Point %s", 
            i, fileName.c_str()), kMagenta+i, dt_cut21_31[0], dt_cut21_31[1], dt_cut21_31[2]);
		auto* h1_dt22_31 = new TH1P(Form("((h1_dt%d_22_31))Delta t [25 ps]@SCI31 - SCI21@TOF Point %s", 
            i, fileName.c_str()), kMagenta+i, dt_cut22_31[0], dt_cut22_31[1], dt_cut22_31[2]);

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
            const auto& sci21 = frs->sci[SCI_21_I];
            const auto& sci22 = frs->sci[SCI_22_I];
			const auto& sci31 = frs->sci[SCI_31_I];
			if(sci22.hits.size() != 1 or sci31.hits.size() != 1 or sci21.hits.size() != 1) 
                continue;
			
			double dt21_31 = sci31.hits[0].t - sci21.hits[0].t;
			double dt22_31 = sci31.hits[0].t - sci22.hits[0].t;
			h1_dt21_31->Fill(dt21_31);
			h1_dt22_31->Fill(dt22_31);

            hist.emplace_back( h1_dt21_31, h1_dt22_31 );

		    auto [res21_31, _ ] = GaussFitMax(*h1_dt21_31, sratio, niter);
		    auto [res22_31, __] = GaussFitMax(*h1_dt22_31, sratio, niter);
            x.push_back( res22_31[1] ); // gauss peak value
        }
    }
    /* TODO! */

    std::time_t now = std::time(nullptr);
    std::tm* tm = std::localtime(&now);
    char buffer[28];
    std::strftime(buffer, sizeof buffer,
        "%Y-%m-%d_%H-%M-%S", tm
    );
	canvas::save_all<canvas::Exe>(save, { buffer });
	WARN("End-of-main");
	rootApp.Run(); return 0;
}


std::istream& operator>>(std::istream& in, Brho& brho) {
    return ::mnd::template operator>> <Brho::SEP>(in, brho.data_);    
}
std::ostream& operator<<(std::ostream& os, const Brho& brho) {
    return os << brho.data_;
}
std::istream& operator>>(std::istream& in, FileBrho& f) {
    char c;
    if(!std::getline(in, f.name, FileBrho::SEP)) {
        WARN("Parsing Brho into string part (file name) failed.\n");
        return in;
    }

    return in >> f.brho;
}
std::ostream& operator<<(std::ostream& os, const FileBrho& f) {
    return os << f.name << ": " << f.brho; 
}
