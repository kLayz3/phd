#include "CMDLineParser.h"

#include "monad/monad.hxx"

#include <regex>
#include <algorithm>

using namespace std;
using CMDLineParser::Mandatory;
const char* Mandatory::def_msg = "Check options by passing --help";

bool CMDLineParser::IsCmdArg(const char* line, int argc, char** argv) {
	char *line1 = (char*)malloc(strlen(line)+3);
	strcpy(line1, "--");
	strcat(line1, line);
	bool retval = 0;
	for(int i(1); i<argc; ++i) {
		if(!strcmp(argv[i], line) || !strcmp(argv[i], line1)) {
			memset(argv[i], '_', strlen(argv[i]));
			retval = 1;
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

template<typename T>
bool 
CMDLineParser::ParseCmdLine(const char* line, T& dest, int argc, char** argv, Mandatory mandatory) {
	using std::string;
	using Vec = std::vector<string>;
	static_assert(std::is_same_v<T, string> || std::is_same_v<T, Vec>,
		"Second argument must be either std::string or std::vector<std::string>>");
	
	/* Regardless if we pass Vector or just String, create a local to reassign later. */
	Vec parsed{};

	cmatch m;
	static const std::regex r(
		R"(^--([^=]+)[=]([^,]+(?:,([^,])+)*)$)"
	);

	bool retval = 0;
	for(int i(1); i<argc; ++i) {
		if(regex_match(argv[i], m, r) and !strcmp(line, m[1].str().c_str())) {
			string list = m[2].str();
			
			size_t pos;
			string split;
			while((pos = list.find(',')) != string::npos) {
				split = list.substr(0, pos);
				parsed.push_back( std::move(split) );
				list.erase(0, pos + 1);
			}
			parsed.push_back( std::move(list) );
			
			// Set argv[i] to be something redundant.
			memset(argv[i], '_', strlen(argv[i]));
			if(retval)
				ERROR("Passed twice argument: \'%s\' in the command line!\n", line);

			retval = 1;
		}
	}
	
	/* Handle the `-tag value` case. */
	for(int i(1); i < argc; ++i) {
		if(argv[i][0] != '-') continue;

		if(!strcmp((char*)(argv[i]+1), line)) {
			memset(argv[i], '_', strlen(argv[i]));
			retval = 1;
			
			/* Eat up all the argv's until we reach another '^-' regex. */
			while((++i) < argc && argv[i][0] != '-') {
				parsed.emplace_back(argv[i]);
				memset(argv[i], '_', strlen(argv[i]));
			}
			--i;
		}
	}

	if(retval) {
		if(parsed.size() == 0)
			ERROR("Saw an argument: \'-%s\', but no followup parameter after.\n", line);

		WARN("Parsed option " EMPH(%s) " with %zu argument%s: ", 
			line, parsed.size(), (parsed.size() > 1) ? "s" : "");
		for(auto& p : parsed) printf(EMPH(%s) " ", p.c_str());
		printf("\n");

		if constexpr(std::is_same_v<T, string>) {
			dest = std::move(parsed.back());
			if(parsed.size() > 1)
				ERROR("(%s): intially parsed as '\%s\', but I found %zu values," 
					" destination marked as single - multiple values not allowed.\n",
					line, parsed[0].c_str(), parsed.size());
		} else { 
			dest = std::move(parsed);
		}
	}
	else if(mandatory.is_it) {
		ERROR("Mandatory argument " EMPH(%s) " not supplied.\n%s%s\n", 
				line, KNRM, mandatory.help_msg ? mandatory.help_msg : mandatory.def_msg);
	}
	return retval;
}

void CMDLineParser::VerifyNoArgumentsLeft(int argc, char** argv) {
	for(int i=1; i < argc; ++i) {
		for(int j=0; j < (int)strlen(argv[i]); ++j) {
			if(argv[i][j] != '_') {
				YELL("Unrecognized or invalid option: " EMPH(%s) "\n", argv[i]);
				printf("Terminating the program.\n");
				exit(11);
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

template bool
CMDLineParser::ParseCmdLine<std::string>(const char*, std::string&, int, char**, Mandatory);

template bool
CMDLineParser::ParseCmdLine<std::vector<std::string>>(const char*, std::vector<std::string>&, int, char**, Mandatory);
