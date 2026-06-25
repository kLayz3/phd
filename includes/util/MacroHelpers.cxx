#include "MacroHelpers.hxx"

std::vector<std::string> ParseFile(const std::string& fileName) {
#ifndef _POSIX_VERSION
#	error "Cannot compile in this function for non- UNIX operating systems!"
#endif
	static constexpr int MAX_BUF_SIZE   = 4096;

	std::string cmd = 
		"gcc -E -P -Werror -undef -x c++ "
		"-fdiagnostics-color=always "
		"-fdiagnostics-show-caret "
		"-ftrack-macro-expansion=0 " 
		"\"" + fileName + "\" " 
		"-o - 2>&1";
	
	std::unique_ptr<FILE, decltype(&pclose)> pipe {popen(cmd.c_str(), "r"), pclose};
	if(!pipe) ERROR("popen failed for file: '\%s\'\n", fileName.c_str());

	std::string text; text.reserve(1024);
	char buf[MAX_BUF_SIZE];

	while(fgets(buf, sizeof(buf), pipe.get())) {
		text += buf;
	}
	
	FILE* pipe_raw = pipe.release();
	int status = pclose(pipe_raw);
	if(WIFSIGNALED(status))
		ERROR("Preprocessor killed by signal %d while parsing '%s'\n",
			WTERMSIG(status), fileName.c_str());
	if(!WIFEXITED(status))
		ERROR("Preproccessor failed:\n " KNRM "%s\n", text.c_str());
	if(WEXITSTATUS(status) != 0)
		ERROR("Preprocessor failed for '%s': " KNRM "\n%s",
			fileName.c_str(), text.c_str());
	
	std::istringstream stream(text);
	std::vector<std::string> lines;

	for(std::string line; std::getline(stream, line); ) {
		if(line.find_first_not_of(" \t\r") != std::string::npos) { // remove empty and whitespace-only lines.
			lines.emplace_back(std::move(line));		
		}
	}
	return lines;
}
