#pragma once

#include <string>
#include <vector>
#include <unordered_set>

namespace CMDLineParser {
	bool IsCmdArg(const char* line, int argc, char** argv);

	struct Mandatory {
		static const char* def_msg;
		bool is_it = false;
		const char* help_msg = nullptr;
		Mandatory() = default;
		Mandatory(bool b) : is_it(b) {}
		static void SetMessage(const char* msg) { def_msg = msg; }
	};
	
	/** 
	 * Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
	/* Or, -tag0 identifier0 -tag1 identifier 1 ... */
	/* If wanting to pass a sequence of identifier to a single tag, the syntax is:
	 *
	 * --tag0=identifier0_0,identifier0_1,identifier0_2 ... (UNIX style) 
	 *  -tag0 identifier0_0 identifier0_1 identifier0_2 ... (Windows style)
	 */
	template<typename T>
	bool ParseCmdLine(const char* line, T& dest, int argc, char** argv, Mandatory mandatory = Mandatory{});

	void VerifyNoArgumentsLeft(int, char**);

	void RemoveCharsFromString(std::string& str, std::string chars); 

	std::vector<std::string> SplitStringToVector(const std::string& str, char delim); 

	std::unordered_set<std::string> SplitStringToSet(const std::string& str, char delim); 
}
