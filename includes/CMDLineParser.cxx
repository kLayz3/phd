#include "CMDLineParser.h"

using namespace std;

bool CMDLineParser::IsCmdArg(const char* line, int argc, char** argv) {
	char *line1 = (char*)malloc(strlen(line)+3);
	strcpy(line1, "--");
	strcat(line1, line);
	for(int i(1); i<argc; ++i) {
		if(!strcmp(argv[i], line) || !strcmp(argv[i], line1)) return 1;
	}
	return 0;
}

/* Cmd args have to be: --tag0=identifier0 --tag1=identifier1 ... */
bool CMDLineParser::ParseCmdLine(const char* line, string& parsed, int argc, char** argv) {
	cmatch m;
	std::regex r("^(--|)([^=]+)[=](.+)$");
	for(int i(1); i<argc; ++i) {
		if(regex_match(argv[i], m, r) && !strcmp(m[2].str().c_str(), line)) {
			parsed = m[3].str(); return 1;
		}
	}   
	return 0;
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

