#include "JSONParser.h"

#include <cstdlib>
#include <stdexcept>
#include <unistd.h>
#include <sys/wait.h>

#ifndef __FILE_NAME__
#	ifdef __unix__
#		define __FILE_NAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#	elif defined(__WIN32) || defined(WIN32)
#		define __FILE_NAME__  (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#	else
#		define __FILE_NAME__ __FILE__
#	endif
#endif

#ifndef KNRM
#	define KNRM "\e[0m"
#endif
#ifndef KCYN
#	define KCYN "\e[0;36m"
#endif
#ifndef KGRN
#	define KGRN "\e[0;32m"
#endif

#ifndef ERROR 
#define ERROR(...) \
do { \
	fprintf(stderr, KGRN "%s" KNRM ":" KCYN "%d" KNRM " => ", __FILE_NAME__, __LINE__); \
	fprintf(stderr, __VA_ARGS__); \
	throw std::runtime_error("JSON Parser: error"); \
} while(0);
#endif

using json = nlohmann::json;
using PipeDeleter = int (*)(FILE*);

json ParseJSON(const std::string& fileName) {
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

	std::unique_ptr<FILE, PipeDeleter> pipe {popen(cmd.c_str(), "r"), pclose};
	if(!pipe) ERROR("popen failed for file: '\%s\'\n", fileName.c_str());

	std::string text;
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

	try {
		return json::parse( text );
	}
	catch(const json::parse_error& e) {
		constexpr std::size_t context = 80;

		// e.byte is 1-based for parse errors.
		const std::size_t pos =
			e.byte > 0 ? e.byte - 1 : 0;

		const std::size_t begin =
			pos > context ? pos - context : 0;

		const std::size_t end =
			std::min(text.size(), pos + context);

		std::cerr
			<< "JSON parse error:\n"
			<< "  " << e.what() << '\n'
			<< "  byte: " << e.byte << '\n'
			<< "  context:\n"
			<< text.substr(begin, end - begin) << '\n';

		throw;
	}
}
