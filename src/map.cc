#include "libs.hh"
#include <algorithm>
#include <csignal>
#include <iostream>
#include "TFile.h"

#include "dbg.hh"
#include "indicators.hh"
#include <csignal>

#include "CMDLineParser.h"
#include "AuxFunctions.hh"
#include "TAnalysisPool.hxx"
#include "TFOOTPedestalProc.h"
#include "TFOOTPedestalCont.h"
#include "TFRSMapProc.h"
#include "TFRSMapCont.h"

using namespace std;
using namespace CMDLineParser;

extern const char* map_help;
void sig_callback_handler(int );

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
	signal(SIGINT | SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;

	srand(time(NULL));

	string pStr, outFile;
	vector<string> fileName{};
	u64 maxEvents = -1;

	CMDLineParser::Mandatory::SetDefMessage(map_help);
	if(IsCmdArg("help", argc, argv)) { cout << def_msg(); return 0; }
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
		int i = ::FindIndex(static_detectors, ID); \
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

	TFRSMapCont frs{};
	frs.Setup(ContainerIO::kOUTPUT, out, h102);

#define MAP_SCI(x, SCI_LABEL) \
	frs.sci[x]._nhit_raw[0] = &sort->tdc_nhit_sc##SCI_LABEL##l; \
	frs.sci[x]._nhit_raw[1] = &sort->tdc_nhit_sc##SCI_LABEL##r; \
	frs.sci[x]._data_raw[0] = &sort->tdc_sc##SCI_LABEL##l[0]; \
	frs.sci[x]._data_raw[1] = &sort->tdc_sc##SCI_LABEL##r[0]; \
	frs.sci[x]._qdc_raw[0]  = &sort->de_##SCI_LABEL##l; \
	frs.sci[x]._qdc_raw[1]  = &sort->de_##SCI_LABEL##r;
	
	MAP_SCI(0, 21);
	MAP_SCI(1, 22);
	MAP_SCI(2, 31);
	MAP_SCI(3, 41);
		
	for(int i=0; i < (int)frs.tpc.size(); ++i) {
		TFRSMapCont::TPC& tpc = frs.tpc[i];
		tpc._tpc_aa = &sort->tpc_a[i][0];
		for(int j=0; j<2; ++j) {
			tpc._tpc_lt[j]  = &sort->tpc_lt[i][j][0];
			tpc._tpc_rt[j]  = &sort->tpc_rt[i][j][0];
			tpc._tpc_ltn[j] = &sort->tpc_nhit_lt[i][j];
			tpc._tpc_rtn[j] = &sort->tpc_nhit_rt[i][j];
		}
		for(int j=0; j<4; ++j) { 
			tpc._tpc_at[j]  = &sort->tpc_dt[i][j][0];
			tpc._tpc_atn[j] = &sort->tpc_nhit_dt[i][j];
		}
	}

	frs.music[0]._music_raw = &sort->music_e1[0];
	frs.music[1]._music_raw = &sort->music_e2[0];
	frs._pattern = &sort->pattern;

	TFOOTPedestalProc::LoadBadStripsFile(PROG_PATH "/params/bad_strips.json");
	auto pool = TAnalysisPool<>()
		.emplace_worker<TFOOTPedestalProc>(foot[0])
		.emplace_worker<TFOOTPedestalProc>(foot[1])
		.emplace_worker<TFOOTPedestalProc>(foot[2])
		.emplace_worker<TFOOTPedestalProc>(foot[3])
		.emplace_worker<TFOOTPedestalProc>(foot[4])
		.emplace_worker<TFOOTPedestalProc>(foot[5])
		.emplace_worker<TFOOTPedestalProc>(foot[6])
		.emplace_worker<TFOOTPedestalProc>(foot[7])
		.emplace_worker<TFRSMapProc>(frs, 0); /* Don't do FRS analysis from first go. */

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
 
	pool.Start();
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar1, ev, nentries);
		pool.AssignWork();
		pool.Await();
	}
	pool.Stop(); bar1.mark_as_completed();
	
	tv.emplace_back(TimePoint("after gped"));
	PrintElapsed<kSECOND>(tv);
	
	/* Perform fitting for the global pedestal calculation. 
	 * Cannot be (obviously) paralellized. */
	for(size_t i=0; i < pool.Size(); ++i) { 
		TFOOTPedestalProc* p = dynamic_cast<TFOOTPedestalProc*>(pool.GetWorker(i));
		if(!p) {
			TFRSMapProc* pfrs = dynamic_cast<TFRSMapProc*>(pool.GetWorker(i));
			if(pfrs) pfrs->do_analysis = 1;
			continue;
		}
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

	pool.Start();
	for(u64 ev = 0; ev < nentries; ++ev) {
		pool.GetEntry(ev);
		PrintProgress(bar2, ev, nentries);
		pool.AssignWork();
		pool.Await();
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
	}
	tv.emplace_back(TimePoint("post fineped fit"));
	PrintElapsed<kSECOND>(tv);

	printf("Total execution time: "); PrintElapsed<kSECOND>(tv.back(), tv.front());

	show_console_cursor(true);

	pool.Write();
	
	out->Close();
	delete h101;
}

const char* map_help =
"\nUsage: ./map <OPT1> <OPT2> ...\n\
Options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input1.root input2.root...   ..Input file(s).\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_cal' suffix.\n\
-help                       ..Print this message to stdout. \n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will go through the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data + perform mapping of the FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";

void sig_callback_handler(int signum) {
	WARN("\nCaught abort/seg signal.\n");
	indicators::show_console_cursor(true);
	exit(signum);
}

