/**
 * Code taken from: https://github.com/jblomer/iotools/blob/master/gen_atlas.cxx
 * This part is omega hacked, because one can't directly write the TFRSXXXEvent object
 * into RNTuple column. Problem is that down the line it inherits from TNamed which is
 * unserializable in modern ROOT parlance.
 */

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TBuffer.h>
#include <TBufferFile.h>

#include "ROOT/RNTuple.hxx"
#include <ROOT/RNTupleModel.hxx>
#include "ROOT/RNTupleWriter.hxx"

#include <cstdio>
#include <filesystem>
#include <getopt.h>
#include <string>

#include <unistd.h>

#include "TString.h"
#include "core/AuxFunctions.hh"

#ifndef GO4_SRC_PATH
#	error "Go4 source path not set."
#endif

#include "TFRSSortEvent.h"

namespace RExp = ROOT::Experimental;

int GetCompressionSettings(std::string shorthand) {
	if (shorthand == "zlib")
		return 101;
	if (shorthand == "lz4")
		return 404;
	if (shorthand == "lzma")
		return 207;
	if (shorthand == "zstd")
		return 505;
	if (shorthand == "none")
		return 0;
	abort();
}

void Usage(const char* progname) {
	printf("Usage: %s args\n"
		"-i <input-file>         .. input file.\n"
		"[-o] <ntuple-path>      .. output destination (dir). Default: $(dirname $input)/rn\n" 
		"[-c] <compression>      .. compression. Default: none\n"
		"\n", progname);
	printf("Will convert a TTree into RNTuple with the same layout and same name.\n");
}

namespace fs = std::filesystem;
using namespace indicators;

int main(int argc, char **argv) {
	std::string inputFile {};
	std::string outputPath = {};
	int compressionSettings = 0;
	std::string compressionShorthand = "none";
	
	std::string treeName = "SortxTree";
	std::string branchName   = "FRSSortEvent.";
	std::string branchNameRN = "FRSSortEvent";

	gSystem->Load(GO4_SRC_PATH "/libGo4UserAnalysis.so");

	int c;
	while ((c = getopt(argc, argv, "hv:i:o:c:t:m")) != -1) {
		switch (c) {
			case 'h':
			case 'v':
				Usage(argv[0]);
				return 0;
			case 'i':
				inputFile = optarg;
				break;
			case 'o':
				outputPath = optarg;
				break;
			case 'c':
				compressionSettings = ::GetCompressionSettings(optarg);
				compressionShorthand = optarg;
				break;
			default:
				fprintf(stderr, "Unknown option: -%c\n", c);
				Usage(argv[0]);
				return 1;
		}
	}

	if(inputFile.empty()) {
		fprintf(stderr, "Input file not provided.\n");
		Usage(argv[0]);
		return 2;
	}
	std::string inputFileName = fs::path{inputFile}.filename();

	if(outputPath.empty()) {
		outputPath = fs::path{inputFile}.parent_path();
		outputPath = (!outputPath.empty()) ? outputPath : ".";

		printf("Output dir not specified. Saving in '%s'\n", outputPath.c_str());
	}

	std::string outputFile = outputPath + "/rn/rn_" + inputFileName;

	unlink(outputFile.c_str());
	printf("Saving into file: %s\n", outputFile.c_str());
	
	TFile fin(inputFile.c_str(), "READ");
	TTree* t = fin.Get<TTree>(treeName.c_str());
	if(!t) {
		fprintf(stderr, "Cannot find the TTree with supplied name.");
		return 3;
	}
	TFRSSortEvent *ev = nullptr;
	t->SetBranchAddress(branchName.c_str(), &ev);
	
	using RawBuf = std::vector<unsigned char>;
	auto model = RExp::RNTupleModel::Create();
	auto fEvt = model->MakeField<RawBuf>(branchNameRN);
	auto writer = RExp::RNTupleWriter::Recreate(std::move(model), treeName, outputFile);

	ProgressBar bar {
		option::BarWidth{50},
		option::Start{"["},
		option::Fill{"-"},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Event-by-event pedestal"},
		option::ForegroundColor{Color::magenta},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	RawBuf buffer;
	u64 nEntries = (u64)t->GetEntries();
	for(Long_t i=0; i<nEntries; ++i) {
		t->GetEntry(i);
		TBufferFile buf(TBuffer::kWrite);
		buf.WriteObjectAny(ev, TFRSSortEvent::Class());
		buffer.assign((unsigned char*)buf.Buffer(),
				(unsigned char*)buf.Buffer() + buf.Length());
		*fEvt = buffer;
		writer->Fill();
		mnd::PrintProgress(bar, (u64)i, nEntries, 500);
	}
	printf("\n");	
	return 0;
}
