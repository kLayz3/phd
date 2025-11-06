/**
 * Code taken from: https://github.com/jblomer/iotools/blob/master/gen_atlas.cxx
 */

#include <ROOT/RNTupleImporter.hxx>

#include <TROOT.h>

#include <cstdio>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <string>

#include <unistd.h>

using RNTupleImporter = ROOT::Experimental::RNTupleImporter;

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
		"[-t] <tree-name>        .. Specify TTree name. Default: 'SortxTree'\n"
		"[-m(t)] <H1 root file>  \n"
		"\n", progname);
	printf("Will convert a TTree into RNTuple with the same layout and same name.\n");
}

namespace fs = std::filesystem;
int main(int argc, char **argv) {
	std::string inputFile {};
	std::string outputPath = {};
	int compressionSettings = 0;
	std::string compressionShorthand = "none";
	
	std::string treeName = "SortxTree";

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
			case 't':
				treeName = optarg;
			case 'm':
				ROOT::EnableImplicitMT();
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
		printf("Output dir not specified. Saving in '%s'\n", outputPath.c_str());
	}

	std::string outputFile = (!outputPath.empty() ? outputPath : "." ) + 
			"/" + "rn" + "_" + inputFileName;

	unlink(outputFile.c_str());
	auto importer = RNTupleImporter::Create(inputFile, treeName, outputFile);
	auto options = importer->GetWriteOptions();
	options.SetCompression(compressionSettings);
	importer->SetWriteOptions(options);
	importer->Import();

	return 0;
}
