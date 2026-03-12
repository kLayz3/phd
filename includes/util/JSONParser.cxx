#include "JSONParser.h"

#include <cstdlib>
#include <stdexcept>
#include <unistd.h>
#include <sys/wait.h>

#ifndef __FILENAME__
#   define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#ifndef KNRM
#define KNRM "\e[0m"
#define KCYN "\e[0;36m"
#define KGRN "\e[0;32m"
#endif

#ifndef ERROR 
#define ERROR(...) \
do { \
	fprintf(stderr, KGRN "%s" KNRM ":" KCYN "%d" KNRM " => ", __FILENAME__, __LINE__); \
	fprintf(stderr, __VA_ARGS__); \
	throw std::runtime_error("JSON Parser: error"); \
} while(0);
#endif

using json = nlohmann::json;

json ParseJSON(const std::string& fileName) {
#ifndef _POSIX_VERSION
#	error "Cannot compile in this function for non- UNIX operating systems!"
#else
	static constexpr int MAX_BUF_SIZE   = 4096;

	std::string cmd = 
		"gcc -E -P -Werror -undef -x c++ "
		"-fdiagnostics-color=always "
		"-fdiagnostics-show-caret "
		"-ftrack-macro-expansion=0 " 
		"\"" + fileName + "\" " 
		"-o - 2>&1";

    FILE* pipe = popen(cmd.c_str(), "r");
    if(!pipe) ERROR("popen failed for file: '\%s\'\n", fileName.c_str());
	
	std::string text;
    char buf[MAX_BUF_SIZE];

    while(fgets(buf, sizeof(buf), pipe)) {
        text += buf;
	}
	
	int status = pclose(pipe);
	if(WIFSIGNALED(status))
		ERROR("Preprocessor killed by signal %d while parsing '%s'\n",
			WTERMSIG(status), fileName.c_str());
	if(!WIFEXITED(status))
		ERROR("Preproccessor failed:\n " KNRM "%s\n", text.c_str());
	if(WEXITSTATUS(status) != 0)
		ERROR("Preprocessor failed for '%s': " KNRM "\n%s",
			fileName.c_str(), text.c_str());

	/* Can throw on bad parse. */
	return json::parse( text );

#endif

}
