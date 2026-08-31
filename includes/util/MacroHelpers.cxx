#include "MacroHelpers.h"

using PipeDeleter = int (*)(FILE*);

template void canvas::save_all<canvas::Macro>(canvas::Extension , std::vector<std::string_view> );
template void canvas::save_all<canvas::Exe  >(canvas::Extension , std::vector<std::string_view> );
template void canvas::save_all<canvas::Macro>(std::vector<canvas::Extension> , std::vector<std::string_view> );
template void canvas::save_all<canvas::Exe  >(std::vector<canvas::Extension> , std::vector<std::string_view> );


void canvas::DumpPrimitives(TVirtualPad* pad, int depth) {
	if(!pad) return;

	TIter next(pad->GetListOfPrimitives());
	while(TObject* obj = next()) {
		std::cout
			<< std::string(depth * 2, ' ')
			<< obj->ClassName()
			<< "  "
			<< obj->GetName()
			<< "  "
			<< obj->GetTitle()
			<< '\n';

		if(auto* subpad = dynamic_cast<TVirtualPad*>(obj))
			DumpPrimitives(subpad, depth + 1);
	}
}

std::vector<TH1*> canvas::GetHistograms(TVirtualPad* pad) {
	std::vector<TH1*> out {};

	TIter next(pad->GetListOfPrimitives());
	while(auto* obj = next()) {
		if(auto* h = dynamic_cast<TH1*>(obj))
			out.push_back(h);

		else if(auto* subpad = dynamic_cast<TVirtualPad*>(obj)) {
			auto sub = GetHistograms(subpad);
			out.insert(out.end(), sub.begin(), sub.end());
		}
	}

	return out;
}

std::vector<std::string> ParseFile(const std::string& fileName) {
#ifndef _POSIX_VERSION
#	error "Cannot compile in this function for non- UNIX operating systems!"
#endif
	static constexpr size_t MAX_BUF_SIZE = (1ULL << 15); // 32 KiB

	std::string cmd = 
		"gcc -E -P -Werror -undef -x c++ "
		"-fdiagnostics-color=always "
		"-fdiagnostics-show-caret "
		"-ftrack-macro-expansion=0 " 
		"\"" + fileName + "\" " 
		"-o - 2>&1";
	
	std::unique_ptr<FILE, PipeDeleter> pipe {popen(cmd.c_str(), "r"), pclose};
	if(!pipe) ERROR("popen failed for file: '\%s\' (%m)\n", fileName.c_str());

	std::string text; text.reserve(1024);
	char buf[MAX_BUF_SIZE];

	while(fgets(buf, sizeof(buf), pipe.get())) {
		text += buf;
	}
	
	FILE* pipe_raw = pipe.release();
	int status = pclose(pipe_raw);
	if(WIFSIGNALED(status))
		ERROR("Preprocessor killed by signal %d while parsing '%s' (%m)\n",
			WTERMSIG(status), fileName.c_str());
	if(!WIFEXITED(status))
		ERROR("Preproccessor failed:\n " KNRM "%s (%m)\n", text.c_str());
	if(WEXITSTATUS(status) != 0)
		ERROR("Preprocessor failed for '%s': " KNRM "%s (%m)\n",
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

std::string ParseFileToString(const std::string& fileName) {
#ifndef _POSIX_VERSION
#	error "Cannot compile in this function for non- UNIX operating systems!"
#endif
	static constexpr size_t MAX_BUF_SIZE = (1ULL << 14);

	std::string cmd = 
		"gcc -E -P -Werror -undef -x c++ "
		"-fdiagnostics-color=always "
		"-fdiagnostics-show-caret "
		"-ftrack-macro-expansion=0 " 
		"\"" + fileName + "\" " 
		"-o - 2>&1";
	
	std::unique_ptr<FILE, PipeDeleter> pipe {popen(cmd.c_str(), "r"), pclose};
	if(!pipe) ERROR("popen failed for file: '\%s\' (%m)\n", fileName.c_str());

	std::string text; text.reserve(1024);
	char buf[MAX_BUF_SIZE];

	while(fgets(buf, sizeof(buf), pipe.get())) {
		text += buf;
	}
	
	FILE* pipe_raw = pipe.release();
	int status = pclose(pipe_raw);
	if(WIFSIGNALED(status))
		ERROR("Preprocessor killed by signal %d while parsing '%s' (%m)\n",
			WTERMSIG(status), fileName.c_str());
	if(!WIFEXITED(status))
		ERROR("Preproccessor failed:\n " KNRM "%s (%m)\n", text.c_str());
	if(WEXITSTATUS(status) != 0)
		ERROR("Preprocessor failed for '%s': " KNRM "%s (%m)\n",
			fileName.c_str(), text.c_str());
	
	return text;
}

using namespace std::literals;

static bool ends_with(std::string_view name, std::string_view extension) {
    if(name.size() < extension.size() ||
       name.substr(name.size() - extension.size()) != extension) 
    {
        return false;
    }
    return true;
}

static bool starts_with(std::string_view name, std::string_view prefix) {
    if(name.size() < prefix.size() ||
       name.substr(0, prefix.size()) != prefix)
    {
        return false;
    }
    return true;
}

std::pair<std::string_view, std::string_view>
mnd::fs::file_number_bounds(const std::string& file) {
    namespace fs = std::filesystem;

    /*
     * Important: get the location of the filename inside the original string,
     * because the returned string_views must refer to `file`, not to some
     * temporary string produced by std::filesystem.
     */
    const fs::path path{file};
    const auto filename = path.filename().string();

    if(filename.size() > file.size())
        throw std::invalid_argument("Invalid file path: " + file);

    const std::size_t filename_pos = file.size() - filename.size();

    std::string_view name {
        file.data() + filename_pos,
        filename.size()
    };

    // Strip ".root".

    if(! ::ends_with(name, FILE_EXTENSION))
        throw std::invalid_argument(
            mnd::msg("Unexpected file extension: expected: \'%s\', received file name: ", FILE_EXTENSION.data(), file.c_str())
        );
    name.remove_suffix(FILE_EXTENSION.size());

    if(! ::starts_with(name, FILE_PREFIX))
        throw std::invalid_argument(
            mnd::msg("Unexpected file prefix, expected: \'%s\', received file name: %s ", FILE_PREFIX.data(), file.c_str())
        );
    name.remove_prefix(FILE_PREFIX.size());

    const auto sep = name.find('_');

    if(sep == std::string_view::npos)
        throw std::invalid_argument(
            mnd::msg("Expected %s_<start>_<end>%s: %s", FILE_PREFIX.data(), FILE_EXTENSION.data(), file.c_str())
        );

    const auto start = name.substr(0, sep);
    const auto end   = name.substr(sep + 1);

    if(start.empty() || end.empty())
        throw std::invalid_argument (
            "Empty file sequence number: " + file
        );

    return {start, end};
}

std::string_view mnd::fs::file_start_number(const std::string& file) {
    return file_number_bounds(file).first;
}

std::string_view mnd::fs::file_end_number(const std::string& file) {
    return file_number_bounds(file).second;
}

std::pair<std::string_view, std::string_view>
mnd::fs::file_number_bounds(const std::vector<std::string>& files) {
    if(files.empty())
        throw std::invalid_argument("Cannot determine bounds of an empty file sequence");

    return {
        file_start_number(files.front()),
        file_end_number(files.back())
    };
}

std::string mnd::fs::file_names_concatenated(const std::vector<std::string>& files) {
    auto bounds = file_number_bounds(files);

    return std::string{FILE_PREFIX} 
        + "_" 
        + std::string{bounds.first} 
        + "_" 
        + std::string{bounds.second};
}
