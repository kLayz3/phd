#include "monad/monad.hxx"

#include <algorithm>
#include <csignal>

#include <csignal>
#include <unistd.h>

#include "CMDLineParser.h"
#include "TFOOTMapProc.h"
#include "TFOOTMapCont.h"
#include "TFRSMapProc.h"
#include "TFRSMapCont.h"

using namespace CMDLineParser;
using namespace mnd;

extern const char* map_help;

/* If this is defined then the original ROOT branch comes from FRS Go4 - Sort step instead of ucesb. */
#if defined(FRS_GO4)
	#pragma message("FRS Go4 already defined somewhere. Careful.")
#else
	#define FRS_GO4
#endif

#define FOOT_ID_0 10
#define FOOT_ID_1 19
#define FOOT_ID_2 17
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

int main(i32 argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;
	CMDLineParser::Mandatory::SetDefMessage(map_help);

	srand(time(NULL));

	std::string pStr, outFile, inFile;
	u64 maxEvents = -1;

	if(IsCmdArg("help", argc, argv)) { std::cout << def_msg(); return 0; }
	
	ParseCmdLine("file", inFile, argc, argv, true /* Necessary arg. */);
	
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = inFile.substr(0, inFile.find('.')) + "_map.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(std::exception const& e) { 
			WARN("Unparsable " EMPH(max-events) " argument to u64. Err: %s\n", e.what()); 
		}
	}

	VerifyNoArgumentsLeft(argc, argv);
	std::vector<TimePoint> tv;

#ifdef FRS_GO4
	TFRSGo4Cont sort{};
#else
	static_assert(false, 
		"Cannot handle non-FRS structures (yet)"); 
#endif

	std::unordered_map<std::string, std::string> info;
	TFOOTMapCont foot[N_FOOT];
#define INIT_FOOT_(ID) \
	{ \
		int i = mnd::FindIndex(static_detectors, ID); \
		if(i < 0) ERROR("Index cannot be found: ID=%d, i=%d", ID, i); \
		TFOOTMapCont& f = foot[i]; \
		info["FOOT_ID"] = #ID; \
		f.Init(info); \
		f.Setup(); \
	}
#define INIT_FOOT(x) INIT_FOOT_(x)

	INIT_FOOT(FOOT_ID_0);
	INIT_FOOT(FOOT_ID_1);
	INIT_FOOT(FOOT_ID_2);
	INIT_FOOT(FOOT_ID_3);
	INIT_FOOT(FOOT_ID_4);
	INIT_FOOT(FOOT_ID_5);
	INIT_FOOT(FOOT_ID_6);
	INIT_FOOT(FOOT_ID_7);

	TFRSMapCont frs{};
	frs.Setup();

	TFOOTMapProc::LoadBadStripsFile(PROG_PATH "/params/bad_strips.json");
	
	u32 n_batch = 25'000;

	auto pool = TAnalysisProcess<>(inFile, outFile, "h102")
		.emplace_process<TFOOTMapProc>( foot[0], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::YES)
		.emplace_process<TFOOTMapProc>( foot[1], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[2], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[3], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[4], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[5], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[6], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFOOTMapProc>( foot[7], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO)
		.emplace_process<TFRSMapProc >( frs,     sort, TFRSMapProc::DoAnalysis::NO)
		.MakePool<1>( 500 );

	/* To register the initial FOOT global pedestals. */
	pool.SendOneBatch( n_batch );
	WARN("Done with batching, going on to process (initial) guess for global pedestal in FOOT\n");

	/* Loop over all workers, get their processes,
	 * and set them correspondingly (FOOT). */
	for(auto& process : pool.GetPool()) {
		auto subprocesses = process.GetProcesses(); 
		/* ^^^^ std::array<TProcessorBase, _> */

		for(TProcessorBase* subproc : subprocesses) { 
			TFOOTMapProc* pfoot = dynamic_cast<TFOOTMapProc*>(subproc);
			if(pfoot) {
				pfoot->CalcGlobalPedestal();
				pfoot->process_type = TFOOTMapProc::kMAIN;
			} else {
				TFRSMapProc* pfrs = dynamic_cast<TFRSMapProc*>(subproc);
				if(pfrs) pfrs->do_analysis = TFRSMapProc::DoAnalysis::YES;
			}
		}
	}
	WARN("Done w/ guessing initial global ped per FOOT, per strip. Proceeding to main mapping step.\n");

	ProgressBar bar {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"-"},
		option::Lead{"@"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Event-by-event pedestal"},
		option::ForegroundColor{Color::blue},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};
	
	tv.emplace_back(TimePoint("start"));

	pool.Start(bar, maxEvents);
	
	pool.Collect();

	/* Perform final fit for the crrected pedestal sigma calculation. */ 
	auto& process = pool.Ref();
	for(TProcessorBase* subproc : process.GetProcesses()) { 
		TFOOTMapProc* pfoot = dynamic_cast<TFOOTMapProc*>(subproc);
		if(pfoot) pfoot->CalcFinalPedestal();
	}
	/* Write gets called here, on pool's dtor. */

	tv.emplace_back(TimePoint("end"));
	PrintElapsed<kSECOND>(tv);

}

const char* map_help =
"\nUsage: ./map <OPT1> <OPT2> ...\n\
Key-value options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input1.root           ..Input file from Go4/UCESB.\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_map' suffix.\n\
-help                       ..Print this message to stdout. \n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will go through the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data + perform mapping of the FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
