#include "TFOOTMapCont.h"
#include "monad/monad.hxx"
#include <csignal>
#include <memory>
#include <string>

#include "TFOOTCalCont.h"
#include "TFRSCalCont.h"
#include "TFRSHitCont.h"
#include "TFRSHitProc.h"
#include "TFOOTHitCont.h"
#include "TFOOTHitProc.h"
#include "util/Verbosity.hxx"
#include "util/CLI.h"

using namespace std::literals;
using namespace mnd;

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	CLI::App app{"This program will do FOOT tracking (and tbd: S2/S3 FRS PID + momentum measurement).\n\
		Always remember: PHYSICS IS FUN <(^.^)>"};

#ifdef MND_FOOTTRACK_DEBUG
    {
        std::string descr = app.get_description();
        descr += "\n" BOLD "Compiled in DEBUG mode. Extra fields appended, and ROOT I/O won't work against "
            "versions with non-debug mode." KNRM;
        app.description(std::move(descr));
    }
#endif
	int verbosity_raw = 0; 
	std::string fileName, outFile, setupFile, footSetupFile;
	u64 maxEvents = -1;
	double kalman_max_cost   = TrackCost::DEFAULT_MAX_CANDIDATE_COST;
	double kalman_max_cost_f = TrackCost::DEFAULT_MAX_FINAL_COST;
	double kalman_cost_cr = TrackCost::DEFAULT_COST_R;
	double kalman_cost_cq = TrackCost::DEFAULT_COST_Q;
	double kalman_cost_ct = TrackCost::DEFAULT_COST_T; 
	bool must_have_upstream_track = false;

	add_logged_option<DisplayDefault::No>(app, "-f,--file", fileName, "Input ROOT file")
		->required()
		->expected(1)
		->check(CLI::ExistingFile)
		->each( [&setupFile, &footSetupFile](const std::string& file_name) {
			auto f = std::make_unique<TFile>(file_name.c_str(), "READ");
			auto* _p1 = f->Get<std::string>("FRS_setup_file");
			if(!_p1 or !_p1->length())
				ERROR("`FRS_setup_file` (std::string) object not found (or is blank) in: %s\n", file_name.c_str());
			setupFile = *_p1;

			_p1 = f->Get<std::string>("FOOT0_setup_file");
			if(!_p1 or !_p1->length())
				ERROR("`FOOT0_setup_file` (std::string) object not found (or is blank) in: %s\n", file_name.c_str());
			footSetupFile = *_p1;
		});


	add_logged_option<DisplayDefault::No>(app, "-o,--output", outFile, 
		"Specify output file name. Default same as the input file with \'_cal\' suffix.")
		->expected(0,1);
	add_logged_option<DisplayDefault::No>(app, "-m,--max-events", maxEvents, 
		"Specify total number of events. Default: all events in the input ROOT file.")
		->check(CLI::PositiveNumber);

	add_logged_option(app, "-u,--upstream", must_have_upstream_track, 
		"Predicate [0=false, 1=true] if upstream track must be present.");
	add_logged_option(app, "-v,--verbose", verbosity_raw, "Verbosity level: 0=silent, 1=info, 2=chatty, 3=spam, 4=infinite")
		->check(CLI::Range(0,4));

	add_logged_option(app, "--cr", kalman_cost_cr, "cr coefficient value.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--cq", kalman_cost_cq, "cq coefficient value.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--ct", kalman_cost_ct, "ct coefficient value.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "-c,--max-cost", kalman_max_cost, "Maximum cost value supplied to Kalman algorithm. Use 0 for infinite cost.")
		->check(CLI::NonNegativeNumber); 
	add_logged_option(app, "-r,--max-cost-final", kalman_max_cost_f, "Maximum cost that the final track can have. Use 0 for infinite cost.")
		->check(CLI::NonNegativeNumber); 

	CLI11_PARSE(app, argc, argv);

	if(kalman_max_cost == 0.0)   kalman_max_cost = HUGE_VAL;
	if(kalman_max_cost_f == 0.0) kalman_max_cost_f = HUGE_VAL;
	if(outFile.empty()) outFile = fileName.substr(0, fileName.find('.')) + "_hit.root"; 
	Verbosity v = *mnd::itov(verbosity_raw);

	srand(time(NULL));
	std::vector<TimePoint> tv;

	TFRSCalCont cfrs;
	cfrs.Setup();

	TFRSHitCont hfrs;
	hfrs.Init( {{"Setup", setupFile }} );
	hfrs.Setup();

	TFOOTCalCont cfoot[N_FOOT_DETECTORS];
	for(int i=0; i<N_FOOT_DETECTORS; ++i) {
		cfoot[i].Init({
			{ "ID", std::to_string(::static_detectors[i] ) }, 
			{ "Setup", footSetupFile }
		});

		cfoot[i].Setup();
	}

	TFOOTHitCont hfoot;
	hfoot.Init( {{"Setup", footSetupFile }} );
	hfoot.Setup();

	/* FRS process must be *before* FOOT process. Must get invoked before FOOT, per event. */
	auto pool = TAnalysisProcess<>(fileName, outFile, "h104")
		.emplace_process<TFRSHitProc>(hfrs    , cfrs, 0xa )
		.emplace_process<TFOOTHitProc>(hfoot,
			cfoot[0], cfoot[1], cfoot[2], cfoot[3], 
			cfoot[4], cfoot[5], cfoot[6], cfoot[7],
			kalman_max_cost, kalman_max_cost_f,
			std::array{ kalman_cost_cr, kalman_cost_cq, kalman_cost_ct },
			must_have_upstream_track, v) 
#ifdef MND_DEBUG_ENABLED 
			.MakePool<1>( 512 );
#elif defined(MND_FOOTTRACK_DEBUG)
			.MakePool<8>( 4092 );
#else
			.MakePool<16>( 4092 );
#endif
	
	ProgressBar bar {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Hit Analysis (per event)"},
		option::ForegroundColor{Color::grey},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector{FontStyle::bold}}
	};

	tv.emplace_back(TimePoint("start"));

	pool.Start(bar, maxEvents);	
	pool.Collect();

	tv.emplace_back(TimePoint("end"));

	PrintElapsed<kSECOND>(std::move(tv));
}
