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
		foot[i].Init( {{"FOOT_ID"s, std::to_string(i)}} );
		foot[i].Setup(ContainerIO::kOUTPUT);
	}
	
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
