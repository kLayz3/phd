#include "monad/monad.hxx"

#include <algorithm>
#include <iostream>
#include <csignal>

#include "util/CMDLineParser.h"
#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "TFOOTCalProc.h"
#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"

using namespace CMDLineParser;
using namespace std::literals;
using namespace mnd;

extern const char* clusterise_help;

int main(int argc, char* argv[]) {
	using namespace indicators;
	signal(SIGINT , sig_callback_handler);
	signal(SIGSEGV, sig_callback_handler);
	auto& def_msg = CMDLineParser::Mandatory::DefMessage;
	CMDLineParser::Mandatory::SetDefMessage(clusterise_help);

	srand(time(NULL));	

	std::string pStr, fileName, outFile;
	u64 maxEvents = -1;

	if(IsCmdArg("help", argc, argv)) { std::cout << def_msg(); return 0; }
	
	ParseCmdLine("file", fileName, argc, argv, true);
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_cal.root"; 
		WARN("No output file specified. Writing to file: %s\n", outFile.c_str());
	}
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
		try { maxEvents = stoi(pStr); }
		catch(std::exception& e) { WARN("Unparsable " EMPH(max-events) " argument to u64"); std::cout << e.what() << std::endl; }
	}
}


const char* clusterise_help =
"\nUsage: ./cal <OPT1> <OPT2> ...\n\
Options can be passed Windows style (-tag value1 value2 ...) or Unix style (--tag=value1,value2,...)\n\
For either single or multiple values.\n\
\n\
-file input.root            ..Input file(s).\n\
-output /PATH/TO/OUT.root   ..Specify output file name. Default same as first input file with '_cal' suffix.\n\
-max-events N               ..Specify how many events to process in the ROOT file. Default all.\n\
-help                       ..Print this message to stdout. \n\
\n\
This program will remake the FOOT clustering on the calibrated ROOT file. FRS data is just forwarded here.\n\
Always remember: PHYSICS IS FUN <(^.^)>\n\n";
