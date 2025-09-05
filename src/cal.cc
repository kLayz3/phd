#include "CMDLineParser.h"
#include "AuxFunctions.hh"
#include "TString.h"
#include "libs.hh"
#include <algorithm>
#include <iostream>
#include "TApplication.h"
#include "TFile.h"

#include "dbg.hh"
#include "indicators.hh"

#include "TAnalysisPool.hxx"
#include "TFOOTPedestalProc.h"
#include "TFOOTPedestalCont.h"

using namespace std;
using namespace CMDLineParser;

extern const char* calibrate_help;

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
	srand(time(NULL));

	string pStr, outFile;
	vector<string> fileName{};
	u64 maxEvents = -1;

	CMDLineParser::Mandatory::SetMessage(calibrate_help);
	if(IsCmdArg("help", argc, argv)) { cout << calibrate_help; return 0; }
	ParseCmdLine("file", fileName, argc, argv, true);
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		auto& ref = fileName[0];
		outFile = ref.substr(0, ref.find('.')) + "_cal.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(exception& e) { WARN("Unparsable " EMPH(max-events) " argument to u64"); cout << e.what() << endl; }
	}

	VerifyNoArgumentsLeft(argc, argv);

	vector<TimePoint> tv;

	TChain* h101 = new TChain(_tree_base_name);
	for(auto& name : fileName) {
		TFile* in = new TFile(name.c_str(), "READ");
		if(!in or in->IsZombie())
			ERROR("Bad input ROOT file: %s\n", name.c_str());
		h101->Add(name.c_str());
	}

	h101->LoadTree(0);

#ifdef FRS_GO4
	TFRSSortEvent* sort{};
	int r = h101->SetBranchAddress(_branch_base_name, &sort);
	if(r != 0) ERROR("SetBranchAddress failed. \'%s\', RC = %d\n", _branch_base_name, r);
#else
	// Pass an h101 generated struct.
	EXT_STR_h101 _sort;
	EXT_STR_h101 *sort = &_sort;
	int r = h101->SetBranchAddress(_branch_base_name, sort);
	if(r != 0) ERROR("SetBranchAddress failed. \'%s\', RC = %d\n", _branch_base_name, r);
#endif

	TFile* out = new TFile(outFile.c_str(), "RECREATE"); 
	TTree* h102 = new TTree("h102", "h102");
	h102->SetAutoFlush(0); h102->SetAutoSave(0);

	std::unordered_map<std::string, std::string> info;
	TFOOTPedestalCont foot[N_FOOT];
#define INIT_FOOT_(ID) \
	{ \
		int i = FindIndex(static_detectors, ID); \
		if(i < 0) ERROR("Index cannot be found: ID=%d, i=%d", ID, i); \
		TFOOTPedestalCont& f = foot[i]; \
		info["FOOT_ID"] = #ID; \
		f.Init(info); \
		f.Setup(ContainerIO::kOUTPUT, out, h102); \
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

	auto pool = TAnalysisPool<>()
		.emplace_worker<TFOOTPedestalProc>(foot[0])
		.emplace_worker<TFOOTPedestalProc>(foot[1])
		.emplace_worker<TFOOTPedestalProc>(foot[2])
		.emplace_worker<TFOOTPedestalProc>(foot[3])
		.emplace_worker<TFOOTPedestalProc>(foot[4])
		.emplace_worker<TFOOTPedestalProc>(foot[5])
		.emplace_worker<TFOOTPedestalProc>(foot[6])
		.emplace_worker<TFOOTPedestalProc>(foot[7]);

	pool.in  = h101; 
	pool.out = h102;

	tv.emplace_back(TimePoint("start"));
	u64 nentries = std::min((u64)pool.GetEntries(), maxEvents);
	dbg("Doing global pedestal analysis...", nentries);
	
	show_console_cursor(false);	
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
 
#ifdef ANALYSIS_MULTITHREADED	
	pool.Start();
#endif
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar1, ev, nentries);
#ifdef ANALYSIS_MULTITHREADED
		pool.AssignWork();
		pool.Await();
#else
		pool.ProcessEntry();
#endif
	}
	pool.Stop(); bar1.mark_as_completed();
	
	tv.emplace_back(TimePoint("after gped"));
	PrintElapsed<kSECOND>(tv);
	
	/* Perform fitting for the global pedestal calculation. 
	 * Cannot be (obviously) paralellized. */
	for(size_t i=0; i < pool.Size(); ++i) { 
		TFOOTPedestalProc* p = dynamic_cast<TFOOTPedestalProc*>(pool.GetWorker(i));
		if(!p) continue;
		p->CalcGlobalPedestal();
		p->process_type = TFOOTPedestalProc::kEPED;
		dbg("Finished with one global pedestal fitting ...", i+1);
	}
	tv.emplace_back(TimePoint("post fit"));
	PrintElapsed<kSECOND>(tv);

	dbg("Doing finer pedestal analysis now...");
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

#ifdef ANALYSIS_MULTITHREADED	
	pool.Start();
#endif
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar2, ev, nentries);
#ifdef ANALYSIS_MULTITHREADED
		pool.AssignWork();
		pool.Await();
#else
		pool.ProcessEntry();
#endif

		pool.FillOutput();
	}
	pool.Stop(); bar2.mark_as_completed();

	tv.emplace_back(TimePoint("after fineped"));
	PrintElapsed<kSECOND>(tv);

	/* Perform final fit for the corrected pedestal sigma calculation. */ 
	for(size_t i=0; i < pool.Size(); ++i) { 
		TFOOTPedestalProc* p = dynamic_cast<TFOOTPedestalProc*>(pool.GetWorker(i));
		if(!p) continue;
		p->CalcFinalPedestal();
		dbg("Finished with one fine pedestal fitting (sigma calc) ...", i+1);
	}
	tv.emplace_back(TimePoint("post fineped fit"));
	PrintElapsed<kSECOND>(tv);

	printf("Total execution time: "); PrintElapsed<kSECOND>(tv.back(), tv.front());

	show_console_cursor(true);

	pool.Write();
	
	out->Close();
	delete h101;
}

const char* calibrate_help =
"\nUsage: ./cal <OPT1> <OPT2> ...\n\
Options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input1.root input2.root...   ..Input file(s).\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_cal' suffix.\n\
-help                       ..Print this message to stdout. \n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will analyse the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
