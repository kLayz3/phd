#include "monad/monad.hxx"

#include <algorithm>
#include <csignal>

#include <csignal>
#include <unistd.h>

#include "TFOOTMapProc.h"
#include "TFOOTMapCont.h"
#include "TFRSMapProc.h"
#include "TFRSMapCont.h"
#include "util/CLI.h"

using namespace mnd;

extern const char* map_help;

/* If this is defined then the original ROOT branch comes from FRS Go4 - Sort step instead of ucesb. */
#if defined(FRS_GO4)
	#pragma message("FRS Go4 already defined somewhere. Careful.")
#else
	#define FRS_GO4
#endif

int main(i32 argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);


	std::string fileName, outFile;
	u64 maxEvents = -1;
	double dt_veto = NAN;

	CLI::App app{"This program will go through the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data "
		"+ perform mapping of the FRS data.\nAlways remember: PHYSICS IS FUN <(^.^)>"};

	add_logged_option<DisplayDefault::No>(app, "-f,--file", fileName, "Input ROOT file from Go4/UCESB")
		->required()
		->expected(1)
		->check(CLI::ExistingFile);
	
	add_logged_option<DisplayDefault::No>(app, "-o,--output", outFile, 
		"Specify output file name. Default same as the input file with \'_map\' suffix.")
		->expected(0,1);

	add_logged_option<DisplayDefault::No>(app, "-m,--max-events", maxEvents, 
		"Specify total number of events. Default: all events in the input ROOT file.")
		->check(CLI::PositiveNumber);

	add_logged_option<DisplayDefault::No>(app, "-d,--foot-dt", dt_veto, 
		"Specify in microseconds FOOT deadtime veto. Consecutive entries with stamp difference <T will be discarded. Default 0")
		->check(CLI::PositiveNumber);

	CLI11_PARSE(app, argc, argv);
	if(outFile.empty()) outFile = fileName.substr(0, fileName.find('.')) + "_map.root"; 
	
	srand(time(NULL));
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
		int i = mnd::FindIndex(::static_detectors, ID); \
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

	TFOOTMapProc::LoadBadStripsFile(PROG_PATH "/params/bad_strips.json");

	TFRSMapCont frs{};
	frs.Setup();

	TTrigMapCont trig{};
	trig.Setup();

	u32 n_batch = 16'384;
	auto pool = TAnalysisProcess<>(fileName, outFile, "h102")
		.emplace_process<TFOOTMapProc>( foot[0], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::YES, dt_veto)
		.emplace_process<TFOOTMapProc>( foot[1], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[2], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[3], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[4], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[5], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[6], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFOOTMapProc>( foot[7], sort, TFOOTMapProc::NBatchPedestal{n_batch}, TFOOTMapProc::CableSwapped::NO , dt_veto)
		.emplace_process<TFRSMapProc >( frs,     sort, TFRSMapProc::DoAnalysis::NO)
		.emplace_process<TTrigMapProc>( trig,    sort, TFRSMapProc::DoAnalysis::NO)
		//.MakePool<4>( n_batch );
		.MakePool<1>( 512 );

	/* To register the initial FOOT global pedestals. */
	pool.SendOneBatch( n_batch );
	WARN("Done with batching, going on to process (initial) guess for global pedestal in FOOT\n");

	/* Loop over all workers, get their processes,
	 * and set them correspondingly (FOOT). FRS+Trig process needs to be switched away from no-op state. */
	for(auto& process : pool.GetPool()) {
		auto subprocesses = process.GetProcesses(); 
		/* ^^^^ std::array<TProcessorBase*, _> */

		for(TProcessorBase* subproc : subprocesses) { 
			TFOOTMapProc* pfoot = dynamic_cast<TFOOTMapProc*>(subproc);
			if(pfoot) {
				pfoot->CalcGlobalPedestal();
				pfoot->process_type = TFOOTMapProc::kMAIN;
			} 
			TFRSMapProc* pfrs = dynamic_cast<TFRSMapProc*>(subproc);
			if(pfrs) pfrs->do_analysis = TFRSMapProc::DoAnalysis::YES;
			
			TTrigMapProc* ptrg = dynamic_cast<TTrigMapProc*>(subproc);
			if(ptrg) ptrg->do_analysis = TFRSMapProc::DoAnalysis::YES;
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
		option::PostfixText{"Event-by-event mapping (FRS) + pedestals (FOOT)"},
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
