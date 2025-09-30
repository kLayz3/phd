#include "libs.hh"
#include <algorithm>
#include <csignal>
#include <iostream>

#include "indicators.hh"
#include <csignal>

#include "CMDLineParser.h"
#include "AuxFunctions.hh"
#include "TAnalysisPool.hxx"
#include "TFOOTMapProc.h"
#include "TFOOTMapCont.h"
#include "TFRSMapProc.h"
#include "TFRSMapCont.h"

using namespace std;
using namespace CMDLineParser;

extern const char* map_help;

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
constexpr i32 N_FOOT = LEN(static_detectors);

int main(i32 argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , util::sig_callback_handler);
	signal(SIGSEGV, util::sig_callback_handler);
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
		outFile = ref.substr(0, ref.find('.')) + "_map.root"; 
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

	std::unordered_map<std::string, std::string> info;
	TFOOTMapCont foot[N_FOOT];
#define INIT_FOOT_(ID) \
	{ \
		int i = ::FindIndex(static_detectors, ID); \
		if(i < 0) ERROR("Index cannot be found: ID=%d, i=%d", ID, i); \
		TFOOTMapCont& f = foot[i]; \
		info["FOOT_ID"] = #ID; \
		f.Init(info); \
		f.Setup(ContainerIO::kOUTPUT, outFile); \
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
	frs.Setup(ContainerIO::kOUTPUT, outFile);

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
		tpc._sci_timerefn = &sort->tpc_nhit_timeref[i];
		tpc._sci_timeref = &sort->tpc_timeref[i][0];
	}

	frs.music[0]._music_raw = &sort->music_e1[0];
	frs.music[1]._music_raw = &sort->music_e2[0];
	frs._pattern = &sort->pattern;

	TFOOTMapProc::LoadBadStripsFile(PROG_PATH "/params/bad_strips.json");
	auto pool = TAnalysisPool<>(h101, outFile, "h102")
		.emplace_worker<TFOOTMapProc>(foot[0])
		.emplace_worker<TFOOTMapProc>(foot[1])
		.emplace_worker<TFOOTMapProc>(foot[2])
		.emplace_worker<TFOOTMapProc>(foot[3])
		.emplace_worker<TFOOTMapProc>(foot[4])
		.emplace_worker<TFOOTMapProc>(foot[5])
		.emplace_worker<TFOOTMapProc>(foot[6])
		.emplace_worker<TFOOTMapProc>(foot[7])
		.emplace_worker<TFRSMapProc>(frs, 0);

	tv.emplace_back(TimePoint("start"));
	u64 nentries = std::min((u64)pool.GetEntries(), maxEvents);
	WARN("Doing global pedestal analysis with (%lu:%s) entries", nentries, util::type_name<u64>().c_str());
	
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
		TFOOTMapProc* p = dynamic_cast<TFOOTMapProc*>(pool.GetWorker(i));
		if(!p) {
			TFRSMapProc* pfrs = dynamic_cast<TFRSMapProc*>(pool.GetWorker(i));
			if(pfrs) pfrs->do_analysis = 1;
			continue;
		}
		p->CalcGlobalPedestal();
		p->process_type = TFOOTMapProc::kEPED;
		WARN("Finished with one global pedestal fitting: %zu\n", i+1);
	}

	tv.emplace_back(TimePoint("post fit"));
	PrintElapsed<kSECOND>(tv);

	WARN("Doing finer pedestal analysis now...\n");
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
		pool.Fill();
	}
	pool.Stop(); bar2.mark_as_completed();

	tv.emplace_back(TimePoint("after fineped"));
	PrintElapsed<kSECOND>(tv);

	/* Perform final fit for the corrected pedestal sigma calculation. */ 
	for(size_t i=0; i < pool.Size(); ++i) { 
		TFOOTMapProc* p = dynamic_cast<TFOOTMapProc*>(pool.GetWorker(i));
		if(!p) continue;
		p->CalcFinalPedestal();
	}
	tv.emplace_back(TimePoint("post fineped fit"));
	PrintElapsed<kSECOND>(tv);

	PrintElapsed<kSECOND>(std::move(tv));

	pool.Write();
	show_console_cursor(true);
}

const char* map_help =
"\nUsage: ./map <OPT1> <OPT2> ...\n\
Key-value options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input1.root input2.root...   ..Input file(s) from Go4/UCESB.\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_map' suffix.\n\
-help                       ..Print this message to stdout. \n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This program will go through the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data + perform mapping of the FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
