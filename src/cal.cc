#include "monad/monad.hxx"

#include <algorithm>
#include <iostream>
#include <csignal>

#include "CMDLineParser.h"
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

#define FOOT_ID_0 10
#define FOOT_ID_1 17
#define FOOT_ID_2 19
#define FOOT_ID_3 20
#define FOOT_ID_4 22
#define FOOT_ID_5 25
#define FOOT_ID_6 23
#define FOOT_ID_7 21

constexpr i32 static_detectors[] = {
	FOOT_ID_0, 
	FOOT_ID_1,
	FOOT_ID_2,
	FOOT_ID_3,
	FOOT_ID_4,
	FOOT_ID_5,
	FOOT_ID_6,
	FOOT_ID_7
};
constexpr i32 N_FOOT = mnd::len(static_detectors);

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;
	CMDLineParser::Mandatory::SetDefMessage(calibrate_help);

	srand(time(NULL));	

	std::string pStr, fileName, outFile;
	u64 maxEvents = -1;

	CMDLineParser::Mandatory::SetDefMessage(calibrate_help);
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

	VerifyNoArgumentsLeft(argc, argv);
	std::vector<TimePoint> tv;

	/* Set up the containers. */
	TFOOTMapCont mfoot[N_FOOT]{}; // input map container.
	for(int i=0; i<N_FOOT; ++i) {
		mfoot[i].Init( {{"FOOT_ID"s, std::to_string(::static_detectors[i])}} );
		mfoot[i].Setup();
	}

	TFRSMapCont mfrs{};
	mfrs.Setup();

	TFRSCalCont cfrs{};
	cfrs.Init( {{"Setup", PROG_PATH "/params/frs_setup.json"}} );
	cfrs.Setup();

	TFOOTCalCont cfoot[N_FOOT]; // output container.
	for(int i=0; i<N_FOOT; ++i) {
		cfoot[i].Init({
			{ "FOOT_ID"s, std::to_string(::static_detectors[i]) }, 
			{ "FOOT_POS"s, std::to_string(i) }
		});
		cfoot[i].Setup();
	}

	/* Set up the process pool. */
	auto pool = TAnalysisProcess<>(fileName, outFile, "h103")
		.emplace_process<TFOOTCalProc>(cfoot[0], mfoot[0], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[1], mfoot[1], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[2], mfoot[2], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[3], mfoot[3], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[4], mfoot[4], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[5], mfoot[5], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[6], mfoot[6], 4, 1)
		.emplace_process<TFOOTCalProc>(cfoot[7], mfoot[7], 4, 1)
		.emplace_process<TFRSCalProc >(cfrs    , mfrs)
		.MakePool<8>( 4096 );
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
-help                       ..Print this message to stdout. \n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will analyse the mapped ROOT file and perform the clustering of the FOOT data + calibrating FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
