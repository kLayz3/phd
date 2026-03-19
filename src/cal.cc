#include "monad/monad.hxx"

#include <algorithm>
#include <iostream>
#include <csignal>

#include "util/CMDLineParser.h"
#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTCalProc.h"
#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"

using namespace CMDLineParser;
using namespace std::literals;
using namespace mnd;

extern const char* calibrate_help;

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;
	CMDLineParser::Mandatory::SetDefMessage(calibrate_help);

	srand(time(NULL));	

	std::string pStr, fileName, outFile, setupFile, footSetupFile;
	u64 maxEvents = -1;

	if(IsCmdArg("help", argc, argv)) { std::cout << def_msg(); return 0; }
	
	ParseCmdLine("file", fileName, argc, argv, true);
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_cal.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(std::exception& e) { WARN("Unparsable " EMPH(max-events) " argument to u64"); std::cout << e.what() << std::endl; }
	}
	if(!ParseCmdLine("setup", setupFile, argc, argv)) {
		setupFile = PROG_PATH "/params/frs_setup.json";
		WARN("FRS setup file not specified. Defaulting to: " EMPH(%s\n), setupFile.c_str());
	}
	if(!ParseCmdLine("foot-setup", footSetupFile, argc, argv)) {
		footSetupFile = PROG_PATH "/params/foot_setup.json";
		WARN("FOOT setup file not specified. Defaulting to: " EMPH(%s\n), footSetupFile.c_str());
	}

	VerifyNoArgumentsLeft(argc, argv);
	std::vector<TimePoint> tv;

	/* Set up the containers. */
	TFOOTMapCont mfoot[N_FOOT]{}; // input map container.
	for(int i=0; i < N_FOOT; ++i) {
		mfoot[i].Init( {{"FOOT_ID"s, std::to_string(::static_detectors[i])}} );
		mfoot[i].Setup();
	}

	TFRSMapCont mfrs{};
	mfrs.Setup();

	TFRSCalCont cfrs{};
	cfrs.Init( {{"Setup", setupFile }} );
	cfrs.Setup();

	TFOOTCalCont cfoot[N_FOOT]; // output container.
	cfoot[0].SetRegisterBox(true);

	for(int i=0; i<N_FOOT; ++i) {
		cfoot[i].Init({
			{ "ID"s, std::to_string( mfoot[i].FOOT_N ) }, 
			{ "Setup"s, footSetupFile }
		});
		cfoot[i].Setup();
	}

	/* Set up the process pool. */
	auto pool = TAnalysisProcess<>(fileName, outFile, "h103")
		.emplace_process<TFOOTCalProc>(cfoot[0], mfoot[0])
		.emplace_process<TFOOTCalProc>(cfoot[1], mfoot[1])
		.emplace_process<TFOOTCalProc>(cfoot[2], mfoot[2])
		.emplace_process<TFOOTCalProc>(cfoot[3], mfoot[3])
		.emplace_process<TFOOTCalProc>(cfoot[4], mfoot[4])
		.emplace_process<TFOOTCalProc>(cfoot[5], mfoot[5])
		.emplace_process<TFOOTCalProc>(cfoot[6], mfoot[6])
		.emplace_process<TFOOTCalProc>(cfoot[7], mfoot[7])
		.emplace_process<TFRSCalProc >(cfrs    , mfrs)
		.MakePool<8>( 4092 );
		//.MakePool<1>( 512 );
		/* Number of subthreads, chunk size. */

	ProgressBar bar {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Clustering & Calibration (per event)"},
		option::ForegroundColor{Color::yellow},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	tv.emplace_back(TimePoint("start"));

	pool.Start(bar, maxEvents);	
	pool.Collect();

	tv.emplace_back(TimePoint("end"));

	PrintElapsed<kSECOND>(std::move(tv));
}

const char* calibrate_help =
"\nUsage: ./cal <OPT1> <OPT2> ...\n\
Options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input.root            ..Input file(s).\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_cal' suffix.\n\
-setup /PATH/TO/FRS.json    ..Specify FRS JSON setup file name.\n\
-foot-setup /PATH/TO/JSON   ..Specify FOOT JSON setup file name.\n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
-help                       ..Print this message to stdout. \n\
\n\
This program will analyse the mapped ROOT file and perform the clustering of the FOOT data + calibrating FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
