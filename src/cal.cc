#include "libs.hh"
#include <algorithm>
#include <iostream>

#include "indicators.hh"
#include <csignal>
#include "CMDLineParser.h"
#include "AuxFunctions.hh"
#include "TAnalysisPool.hxx"
#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTCalProc.h"
#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"

using namespace CMDLineParser;
using namespace std::literals;

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
constexpr i32 N_FOOT = LEN(static_detectors);

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , util::sig_callback_handler);
	signal(SIGSEGV, util::sig_callback_handler);
	CMDLineParser::Mandatory::SetDefMessage(calibrate_help);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;

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
	TFOOTMapCont foot[N_FOOT]{}; // input container.
	for(int i=0; i<N_FOOT; ++i) {
		foot[i].Init( {{"FOOT_ID"s, std::to_string(::static_detectors[i])}} );
		foot[i].Setup(ContainerIO::kINPUT_FULL, fileName);
	}

	TFRSMapCont frs{};
	frs.Setup(ContainerIO::kINPUT_RNONLY, fileName);

	TFRSCalCont cfrs{};
	cfrs.Setup(ContainerIO::kOUTPUT, outFile); 
	cfrs.Init( {{"Setup", PROG_PATH "/params/frs_setup.json"}} );

	TFOOTCalCont cfoot[N_FOOT]; // output container.
	for(int i=0; i<N_FOOT; ++i) {
		cfoot[i].Init({
			{ "FOOT_ID"s, std::to_string(::static_detectors[i]) }, 
			{ "FOOT_POS"s, std::to_string(i) }
		});
		cfoot[i].Setup(ContainerIO::kOUTPUT, outFile);
	}

	/* Set up the processing pool. */
	auto pool = TAnalysisPool<>()
		.emplace_worker<TFOOTCalProc>(foot[0], cfoot[0], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[1], cfoot[1], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[2], cfoot[2], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[3], cfoot[3], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[4], cfoot[4], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[5], cfoot[5], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[6], cfoot[6], 4, 1)
		.emplace_worker<TFOOTCalProc>(foot[7], cfoot[7], 4, 1)
		.emplace_worker<TFRSCalProc>(frs, cfrs);
	pool.SetInput(fileName);
	pool.SetOutput(outFile, "h103");

	tv.emplace_back(TimePoint("start"));
	u64 nentries = std::min((u64)pool.GetEntries(), maxEvents);

	ProgressBar bar {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Clustering (per event)"},
		option::ForegroundColor{Color::yellow},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	show_console_cursor(false);
	pool.Start();
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar, ev, nentries);
		pool.AssignWork();
		pool.Await();
		pool.Fill();
	}
	pool.Stop(); bar.mark_as_completed();
	tv.emplace_back(TimePoint("end"));

	PrintElapsed<kSECOND>(std::move(tv));

	show_console_cursor(true);
	pool.Write();
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
