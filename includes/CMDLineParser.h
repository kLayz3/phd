#pragma once

#include <cstdio>
#include <regex>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>

using namespace std;

namespace CMDLineParser {
	bool IsCmdArg(const char* line, int argc, char** argv);

	/* Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
	bool ParseCmdLine(const char* line, std::string& parsed, int argc, char** argv); 

	void RemoveCharsFromString(std::string& str, std::string chars); 

	vector<string> SplitStringToVector(const std::string& str, char delim); 

	unordered_set<string> SplitStringToSet(const std::string& str, char delim); 
}
