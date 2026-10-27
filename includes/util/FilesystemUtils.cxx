#include "FilesystemUtils.h"

#include <cerrno>
#include <cctype>
#include <cstdio>
#include <cstring>

#include "monad/accessories/typedefs.hxx"

using PipeDeleter = int (*)(FILE*);

std::string mnd::fs::parse_file_to_string(const std::string& fileName) {
#ifndef _POSIX_VERSION
#	error "Cannot compile in this function for non- UNIX operating systems!"
#endif
	
	if(fileName.empty() || fileName.find('\0') != std::string::npos)
		throw std::invalid_argument("Invalid config filename");

	// Absolute paths cannot be interpreted as GCC command-line options.
	const auto path = std::filesystem::absolute(fileName).string();
	std::string quoted = "'";
	for(char c : path) {
		if(c == '\'') quoted += "'\\''";
		else          quoted += c;
	}

	quoted += '\'';
	const std::string cmd =
		"gcc -E -P -Werror -undef -x c++ "
		"-fdiagnostics-color=always "
		"-fdiagnostics-show-caret "
		"-ftrack-macro-expansion=0 "
		 + quoted + " -o - 2>&1";
	
	std::unique_ptr<FILE, PipeDeleter> pipe {popen(cmd.c_str(), "r"), pclose};
	if(!pipe) ERROR("popen failed for file: '%s' (%m)\n", fileName.c_str());

	std::string text; text.reserve(1024);
	char buf[1U << 14];

	while(fgets(buf, sizeof(buf), pipe.get())) {
		text += buf;
	}
	
	int status = pclose(pipe.release());
	if(status == -1)
		MND_THROW("pclose failed for '%s': %s\n", fileName.c_str(), std::strerror(errno));
	if(WIFSIGNALED(status))
		MND_THROW("Preprocessor killed by signal %d while parsing '%s' (%m)\n",
			WTERMSIG(status), fileName.c_str());
	if(!WIFEXITED(status))
		MND_THROW("Preproccessor failed:\n " KNRM "%s (%m)\n", text.c_str());
	if(WEXITSTATUS(status) != 0)
		MND_THROW("Preprocessor failed for '%s': " KNRM "%s (%m)\n",
			fileName.c_str(), text.c_str());
	
	return text;
}

std::vector<std::string> mnd::fs::parse_file_to_lines(const std::string& fileName) {
	auto text = parse_file_to_string(fileName);
	
	std::istringstream stream(text);
	std::vector<std::string> lines;

	for(std::string line; std::getline(stream, line); ) {
		if(line.find_first_not_of(" \t\r") != std::string::npos) { // remove empty and whitespace-only lines.
			lines.emplace_back(std::move(line));
		}
	}
	return lines;
}

std::filesystem::path mnd::fs::resolve_maybe_symlink(const std::filesystem::path& path) {
	if(!std::filesystem::is_symlink(path))
		return path;

	auto target = std::filesystem::read_symlink(path);
	return target.is_absolute() ? target : path.parent_path() / target;
}

