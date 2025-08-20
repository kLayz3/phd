#include "CMDLineParser.h"
#include <cstdio>
#include <cstring>
#include <regex>
#include <algorithm>
#include "libs.hh"

using namespace std;

bool CMDLineParser::IsCmdArg(const char* line, int argc, char** argv) {
	char *line1 = (char*)malloc(strlen(line)+3);
	strcpy(line1, "--");
	strcat(line1, line);
	bool retval = 0;
	for(int i(1); i<argc; ++i) {
		if(!strcmp(argv[i], line) || !strcmp(argv[i], line1)) {
			memset(argv[i], '_', strlen(argv[i]));
			retval = 1;
			WARN("Parsed " EMPH(%s) "\n", line);
		}
	}

	for(int i(1); i<argc; ++i) {
		if(argv[i][0] != '-') continue;
		if(!strcmp((char*)(argv[i]+1), line)) {
			/* Next (possible) arg must be a positional arg starting with '-', and not be a value. */
			if(i == argc-1 or argv[i+1][0] == '-') {
				memset(argv[i], '_', strlen(argv[i]));
				retval = 1;
			}
		}
	}
	if(retval) WARN("Parsed " EMPH(%s) "\n", line);
	free(line1); return retval;
}

/* Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
bool CMDLineParser::ParseCmdLine(const char* line, string& parsed, int argc, char** argv) {
	cmatch m;
	std::regex r("^(?:--)([^=]+)[=](.+)$");
	bool retval = 0;
	for(int i(1); i<argc; ++i) {
		if(regex_match(argv[i], m, r) && !strcmp(m[1].str().c_str(), line)) {
			parsed = m[2].str();
			// Set argv[i] to be something redundant.
			memset(argv[i], '_', strlen(argv[i]));
			retval = 1;
		}
	}
	
	/* Handle the `-tag value` case. */
	for(int i(1); i<argc-1; ++i) {
		if(argv[i][0] != '-') continue;

		if(!strcmp((char*)(argv[i]+1), line) and argv[i+1][0] != '-') {
			parsed = std::string(argv[i+1]);

			// Set argv[i] and argv[i+1] to be something redundant.
			memset(argv[i], '_', strlen(argv[i]));
			memset(argv[i+1], '_', strlen(argv[i+1]));
			retval = 1;
		}
	}
	if(retval) WARN("From " EMPH(%s) " parsed " EMPH(%s) "\n", line, parsed.c_str());
	  return retval;
}

void CMDLineParser::VerifyNoArgumentsLeft(int argc, char** argv) {
	for(int i=2; i < argc; ++i) {
		for(int j=0; j < (int)strlen(argv[i]); ++j) {
			if(argv[i][j] != '_') {
				YELL("Unrecognized or invalid option: " EMPH(%s) "\n", argv[i]);
				printf("Terminating the program.\n");
				exit(1111);
			}
		}
	}
}

void CMDLineParser::RemoveCharsFromString(string& str, string chars) {
	for(char c : chars) {
		str.erase(std::remove(str.begin(), str.end(), c), str.end());
	}
}

vector<string> CMDLineParser::SplitStringToVector(const string& str, char delim) {
	std::stringstream iss(str);
	vector<string> parts;
	string part;
	while(std::getline(iss, part, delim)) {
		parts.push_back(part);
	}
	return parts;
}

unordered_set<string> CMDLineParser::SplitStringToSet(const string& str, char delim) {
	std::stringstream iss(str);
	unordered_set<string> parts;
	string part;
	while(std::getline(iss, part, delim)) {
		parts.insert(part);
	}
	return parts;
}

