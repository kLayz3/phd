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

#if 1
#define ANALYSIS_MULTITHREADED
#else
#warning "Running single-threaded. Possibly slower."
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
	
	string pStr, fileName, outFile;
	u64 maxEvents = -1;

	if(argc < 2) {	
		YELL("Must supply a file argument!\n");
		printf("%s", clusterize_help);
		return 0;
	}
}

const char* clusterize_help =
"\nUsage: ./cal <OPT1> <OPT2> ...\n\
\n\
[--file=]inputName.root      ..Input file.\n\
--output=/PATH/TO/OUT.root   ..Specify output file name. Default same as input file with '_subtr' suffix.\n\
--help                       ..Print this message to stdout. \n\
--max-events=N               ..Specify how many events to process in the ROOT file. Default all.\n\
\n\
This script will analyse the raw (sorted) ROOT file and do the full pedestal analysis of the FOOT data.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
