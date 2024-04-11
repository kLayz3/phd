#include "CMDLineParser.h"
#include "AuxFunctions.h"
#include "../libs.hh"

using namespace std;
using namespace CMDLineParser;

const char* unpack_event = "FRSUnpackEvent";

void SetAllBranchAddress(TTree* h101, uint32_t* FOOT, uint32_t* FOOTE[640]) {
	if(!h101 || h101->IsZombie()) return;
	for(int i=0; i<30; ++i) {
		h101->SetBranchAddress(TString::Format("%s.FOOT%d", unpack_event, i), &FOOT[i]);
		h101->SetBranchAddress(TString::Format("%s.FOOT%dE", unpack_event, i), &FOOTE[i]);
	}
}
	
void SubtractPed(const char* fileName, const char* outFile, u64 firstEvent=0, u64 maxEvents=0) {	
	TFile* in = new TFile(fileName, "READ");
	if(!in || in->IsZombie()) {cerr << "Can't open rootfile with name: " << fileName << "\n"; exit(EXIT_FAILURE);}
	
	TTree* h101 = static_cast<TTree*>(in->Get("UnpackxTree"));
	
	/* Read Containers */
	uint32_t FOOT[30];
	uint32_t FOOTE[30][640];

	/* Write containers */
	double pedestal[30][640] = {0};
	double sigPedestal[30][640] = {0};
}

auto main(int argc, char* argv[]) -> int {
    u64 firstEvent = 0;
	u64 maxEvents = 0;
	if(!ParseCmdLine("file", fileName, argc, argv)) {
		cerr << "No file specified!\n";
		return 0;
	}
	if(!ParseCmdLine("output", outFile, argc, argv)) {
		outFile = fileName.substr(0, fileName.find('.')) + "_subtr.root"; 
		cout << "No output file specified. Writing into file: " << outFile << endl;
	}
	
	if(ParseCmdLine("first-event", pStr, argc, argv)) {
        try {
            firstEvent = stoul(pStr);
            printf("Starting from ev#: %ld\n", firstEvent);
        }
        catch(exception& e) {}
    }
	
	if(ParseCmdLine("max-events", pStr, argc, argv)) {
        try {
            maxEvents = stoul(pStr);
            printf("Max events: %ld\n", maxEvents);
        }
        catch(exception& e) {}
    }

	
}

