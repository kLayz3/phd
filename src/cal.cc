#include "CMDLineParser.h"
#include "AuxFunctions.hh"
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

using namespace std;
using namespace CMDLineParser;

extern const char* calibrate_help;

#if 1
#define ANALYSIS_MULTITHREADED
#else
#warning "Running single-threaded. Possibly slower."
#endif


/* If this is defined then the original ROOT branch comes from FRS Go4 - Sort step instead of ucesb. */
#if defined(FRS_GO4)
	#pragma message("FRS Go4 already defined somewhere. Careful.")
#else
	#define FRS_GO4
#endif

#include "TFRSSortEvent.h"

constexpr const char* _tree_base_name = 
#ifdef FRS_GO4
	 "SortxTree"
#else
	"h101"
#endif 
	;
constexpr const char* _branch_base_name =
#ifdef FRS_GO4
	 "FRSSortEvent."
#else
	""
#endif 
	;

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

auto main(i32 argc, char* argv[]) -> i32 {
	using namespace indicators;
	show_console_cursor(false);	

	string pStr, fileName, outFile;
	u64 maxEvents = -1;

	if(argc < 2) {	
		YELL("Must supply a file argument!\n");
		printf("%s", calibrate_help);
		return 0;
	}

	if(IsCmdArg("help", argc, argv)) { cout << calibrate_help; return 0; }
	
	if(!ParseCmdLine("file", fileName, argc, argv)) {
		fileName = std::string(argv[1]);
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_cal.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(exception& e) { WARN("Unparsable " EMPH(max-events) " argument to u64"); cout << e.what() << endl; }
	}

	VerifyNoArgumentsLeft(argc, argv);

	TApplication* app = new TApplication("myApp", 0, 0);
	vector<TimePoint> tv;

	TFile* in = new TFile(fileName.c_str(), "READ");
	if(!in or in->IsZombie()) {
		WARN("Bad input ROOT file: %s\n", fileName.c_str());
		exit(-2);
	}
	TTree* h101 = dynamic_cast<TTree*>(in->Get(_tree_base_name));
	if(!h101 or h101->IsZombie()) {
		YELL("TTree static_cast is somehow nullptr?\n");
		exit(-3);
	}

	TFRSSortEvent* sort;
	int r = h101->SetBranchAddress(_branch_base_name, &sort);
	if(r != 0) ERROR("SetBranchAddress failed. \'%s\', RC = %d\n", _branch_base_name, r);

	TFile* out = new TFile(outFile.c_str(), "RECREATE"); 
	TTree* h102 = new TTree("h102", "h102");

	std::unordered_map<std::string, std::string> info;
	TFOOTPedestalCont foot[N_FOOT];
#define INIT_FOOT_(ID) \
	{ \
		int i = FindIndex(static_detectors, ID); \
		if(i < 0) ERROR("Index cannot be found: ID=%d, i=%d", ID, i); \
		TFOOTPedestalCont& f = foot[i]; \
		info["FOOT_ID"] = #ID; \
		f.Init(info); \
		h102->Branch(foot[i].GetName(), &foot[i]); \
		f._FOOT = &sort->FOOT##ID; \
		f._FOOTE = sort->FOOT##ID##E; \
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

	TAnalysisPool<8> pool;
	pool.out = h102;
	pool.in = h101; 
	FOR(i, N_FOOT) pool.AddOwnedWorker(new TFOOTPedestalProc(foot[i]));

	pool.FinalizeInit();
	tv.emplace_back(TimePoint("start"));
	dbg("Doing global pedestal analysis...");
	u64 nentries = std::min((u64)pool.GetEntries(), maxEvents);
	
	ProgressBar bar1 {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Global Pedestal"},
		option::ForegroundColor{Color::green},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	pool.Start();
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar1, ev, nentries);
#ifdef ANALYSIS_MULTITHREADED
		pool.AssignWork();
		pool.Await();
#else
		for(int i=0; i < (int)pool.Size(); ++i)
			pool.GetWorker(i)->ProcessEntry();
#endif
	}
	pool.Stop(); bar1.mark_as_completed();
	
	tv.emplace_back(TimePoint("after gped"));
	PrintElapsed<kSECOND>(tv);

	/* Set the worker to calculate entry-specific pedestal.
	 * Also save the global pedestal calculation. */
	for(size_t i=0; i < pool.Size(); ++i) { 
		TFOOTPedestalProc* p = dynamic_cast<TFOOTPedestalProc*>(pool.GetWorker(i));
		if(!p) continue;
		p->CalcGlobalPedestal();
		p->do_global_pedestal = false;
	}

	dbg("Doing finer pedestal analysis...");
	ProgressBar bar2 {
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

	pool.Start();
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar2, ev, nentries);
#ifdef ANALYSIS_MULTITHREADED
		pool.AssignWork();
		pool.Await();
#else
		for(int i=0; i < (int)pool.Size(); ++i)
			pool.GetWorker(i)->ProcessEntry();
#endif

		pool.FillOutput();
	}
	pool.Stop(); bar2.mark_as_completed();

	tv.emplace_back(TimePoint("after fineped"));
	PrintElapsed<kSECOND>(tv);

	printf("Total execution time: "); PrintElapsed<kSECOND>(tv.back(), tv.front());

	show_console_cursor(true);

	out->Write();
	out->Close();

	app->Run();
	
	in->Close();
}

const char* calibrate_help =
"\nUsage: ./clusterise <OPT1> <OPT2> ...\n\
\n\
[--file=]inputName.root      ..Input file.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with '_subtr' suffix.\n\
--help                       ..Print this message to stdout. \n\
--max-events=N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
