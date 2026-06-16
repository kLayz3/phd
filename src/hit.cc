#include "TFOOTMapCont.h"
#include "monad/monad.hxx"
#include <csignal>
#include <exception>
#include <memory>
#include <string>

#include "util/CMDLineParser.h"
#include "TFOOTCalCont.h"
#include "TFOOTCalProc.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"
#include "TFRSHitCont.h"
#include "TFRSHitProc.h"
#include "TFOOTHitCont.h"
#include "TFOOTHitProc.h"
#include "util/Verbosity.hxx"

using namespace CMDLineParser;
using namespace std::literals;
using namespace mnd;

extern const char* hit_help;

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;
	CMDLineParser::Mandatory::SetDefMessage(hit_help);

	if(IsCmdArg("help", argc, argv)) { std::cout << def_msg(); return 0; }
	
	std::string pStr, fileName, outFile, setupFile, footSetupFile;
	
	u64 maxEvents = -1;
	double foot_qt = TFOOTHitProc::DEFAULT_MAX_Q_TOLERANCE;
	double foot_mc = TFOOTHitProc::DEFAULT_MAX_COST;
	Verbosity v = Verbosity::SILENT;

	ParseCmdLine("file", fileName, argc, argv, true);
	{
		auto f = std::make_unique<TFile>(fileName.c_str(), "READ");
		auto* _p1 = f->Get<std::string>("FRS_setup_file");
		if(!_p1 or !_p1->length())
			ERROR("`FRS_setup_file` (std::string) object not found (or is blank) in: %s\n", fileName.c_str());
		setupFile = *_p1;

		_p1 = f->Get<std::string>("FOOT0_setup_file");
		if(!_p1 or !_p1->length())
			ERROR("`FOOT0_setup_file` (std::string) object not found (or is blank) in: %s\n", fileName.c_str());
		footSetupFile = *_p1;
	}
	
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_hit.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try {
			maxEvents = std::stoull(pStr);
		} catch(const std::exception& e) {
			ERROR("Conversion \'%s\' to unsigned long long failed. Reason: %s", pStr.c_str(), e.what());
		}
	}
	if(ParseCmdLine("v", pStr, argc, argv)) {
		int val;
		try {
			val = std::stoi(pStr);
		} catch(std::exception& e) {
			ERROR("Value parsed from `v`: \'%s\' not convertible to int.\n", pStr.c_str());
		}
		auto ov = mnd::itov(val);
		v = ov ? *ov : Verbosity::CHATTY;
	}

	VerifyNoArgumentsLeft(argc, argv);
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
			foot_qt, foot_mc, v) 
		.MakePool<4>( 4092 );
		//.MakePool<1>( 512 );
	
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
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	tv.emplace_back(TimePoint("start"));

	pool.Start(bar, maxEvents);	
	pool.Collect();

	tv.emplace_back(TimePoint("end"));

	PrintElapsed<kSECOND>(std::move(tv));
}


const char* hit_help =
"\nUsage: ./hit <OPT1> <OPT2> ...\n\
Options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input.root            ..Input file(s).\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_cal' suffix.\n\
-max-events N               ..Specify total number of events. Default: all events in the input ROOT file.\n\
-v    [0,1,2]               ..Specify verbosity. Default 0 (silent).\n\
-help                       ..Print this message to stdout. \n\
\n\
This program will do FOOT tracking (and tbd: S2/S3 FRS PID + momentum measurement).\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
