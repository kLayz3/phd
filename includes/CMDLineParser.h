#pragma once

#include <string>
#include <vector>
#include <unordered_set>

namespace CMDLineParser {
	bool IsCmdArg(const char* line, int argc, char** argv);

	/* Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
	bool ParseCmdLine(const char* line, std::string& parsed, int argc, char** argv); 
	void VerifyNoArgumentsLeft(int, char**);
	void RemoveCharsFromString(std::string& str, std::string chars); 

	std::vector<std::string> SplitStringToVector(const std::string& str, char delim); 

	std::unordered_set<std::string> SplitStringToSet(const std::string& str, char delim); 
}
