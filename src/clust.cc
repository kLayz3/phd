#include "CMDLineParser.h"
#include "AuxFunctions.hh"
#include "TContainer.h"
#include "TString.h"
#include "libs.hh"
#include <algorithm>
#include <iostream>
#include "TApplication.h"
#include "TFile.h"
#include "TSystem.h"
#include "dbg.hh"

#include "indicators.hh"

#include "TAnalysisPool.hxx"
#include "TFOOTPedestalProc.h"
#include "TFOOTPedestalCont.h"
#include "TFOOTCalProc.h"
#include "TFOOTCalCont.h"

using namespace CMDLineParser;
using namespace std::literals;

#if !defined(ANALYSIS_MULTITHREADED)
	/* Default build: enable multithread. */
	#if 1
		#define ANALYSIS_MULTITHREADED
	#else
		#warning "Running single-threaded. Possibly slower for complex `ProcessEntry` calls!"
	#endif
#endif

#if defined(ANALYSIS_SINGLETHREADED)
	#undef ANALYSIS_MULTITHREADED
	#warning "Running single-threaded. Possibly slower for complex `ProcessEntry` calls!"
#endif

extern const char* clusterize_help;

#define FOOT_ID_0 25
#define FOOT_ID_1 23
#define FOOT_ID_2 22
#define FOOT_ID_3 21
#define FOOT_ID_4 20
#define FOOT_ID_5 19
#define FOOT_ID_6 17
#define FOOT_ID_7 10

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

auto main(int argc, char* argv[]) -> i32 {
	using namespace indicators;
	show_console_cursor(false);
	srand(time(NULL));	

	std::string pStr, fileName, outFile;
	u64 maxEvents = -1;

	if(argc < 2) {	
		YELL("Must supply a file argument!\n");
		printf("%s", clusterize_help);
		return 0;
	}
	if(IsCmdArg("help", argc, argv)) { std::cout << clusterize_help; return 0; }
	
	if(!ParseCmdLine("file", fileName, argc, argv)) {
		fileName = std::string(argv[1]);
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_cal.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(std::exception& e) { WARN("Unparsable " EMPH(max-events) " argument to u64"); std::cout << e.what() << std::endl; }
	}

	std::vector<TimePoint> tv;

	VerifyNoArgumentsLeft(argc, argv);
	
	TFile* in = new TFile(fileName.c_str(), "READ");
	if(!in or in->IsZombie())
		ERROR("Bad input ROOT file: %s\n", fileName.c_str());
	TTree* h102 = dynamic_cast<TTree*>(in->Get("h102"));
	if(!h102 or h102->IsZombie())
		ERROR("TTree cast is somehow nullptr?\n");

	TFile* out = new TFile(outFile.c_str(), "RECREATE"); 
	TTree* h103 = new TTree("h103", "h103");
	h103->SetAutoFlush(0); h103->SetAutoSave(0);

	std::unordered_map<std::string, std::string> info;
	TFOOTPedestalCont foot[N_FOOT]; // input container.
	for(int i=0; i<N_FOOT; ++i) {
		foot[i].Init( {{"FOOT_ID"s, std::to_string(::static_detectors[i])} } );
		foot[i].Setup(ContainerIO::kINPUT);
	}
	
	TFOOTCalCont cfoot[N_FOOT]; // output container.
	for(int i=0; i<N_FOOT; ++i) {
		cfoot[i].Init({
			{ "FOOT_ID"s, std::to_string(::static_detectors[i]) }, 
			{ "FOOT_POS"s, std::to_string(i) }
		});
		cfoot[i].Setup(ContainerIO::kOUTPUT);
	}
	auto pool = TAnalysisPool<>()
		.emplace_worker<TFOOTCalProc>(foot[0], cfoot[0], 6, 2)
		.emplace_worker<TFOOTCalProc>(foot[1], cfoot[1], 6, 1)
		.emplace_worker<TFOOTCalProc>(foot[2], cfoot[2], 6, 3)
		.emplace_worker<TFOOTCalProc>(foot[3], cfoot[3], 6, 1.5)
		.emplace_worker<TFOOTCalProc>(foot[4], cfoot[4], 5, 2)
		.emplace_worker<TFOOTCalProc>(foot[5], cfoot[5], 5, 1)
		.emplace_worker<TFOOTCalProc>(foot[6], cfoot[6], 5, 3)
		.emplace_worker<TFOOTCalProc>(foot[7], cfoot[7], 5, 1.5);
	
	pool.in = h102;
	pool.out = h103;

	tv.emplace_back(TimePoint("start"));
	dbg("Doing the cluster analysis...");
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

#ifdef ANALYSIS_MULTITHREADED	
	pool.Start();
#endif
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar, ev, nentries);
#ifdef ANALYSIS_MULTITHREADED
		pool.AssignWork();
		pool.Await();
#else
		pool.ProcessEntry();
#endif
		pool.FillOutput();
	}
	pool.Stop(); bar.mark_as_completed();
	tv.emplace_back(TimePoint("end"));

	printf("Total execution time: "); PrintElapsed<kSECOND>(tv.back(), tv.front());

	show_console_cursor(true);
	pool.Write();

	out->Close();
	in->Close();
}

const char* clusterize_help =
"\nUsage: ./cal <OPT1> <OPT2> ...\n\
\n\
[--file=]inputName.root      ..Input file.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with '_subtr' suffix.\n\
--help                       ..Print this message to stdout. \n\
--max-events=N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will analyse the calibrated ROOT file and perform the clustering of the FOOT data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
