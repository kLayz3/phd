#include "MacroHelpers.h"
#include "FromChars.h"
#include <string_view>

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

std::vector<std::string> mnd::ParseFile(const std::string& fileName) {
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

std::string mnd::ParseFileToString(const std::string& fileName) {
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

static const std::regex re {mnd::fs::filename_pattern};

std::tuple<std::string, u32, u32>
mnd::fs::file_info(std::string_view file) {
    namespace fs = std::filesystem;

    /* Important: get the location of the filename inside the original string,
     * because the returned string_view must refer to `file`, not to some
     * temporary string produced by std::filesystem. */
    const fs::path path{file};
    const std::string name = path.filename().string();

    if(name.size() > file.size()) {
        MND_THROW("Invalid file path? %*s", (int)file.size(), file.data());
	}

	std::smatch match;
    if(!std::regex_match(name, match, re)) {
		MND_THROW("mnd::fs::file_info: provided file name %s%s%s does not "
			"match the regular expression: %s%s%s",
			BOLD, name.c_str(), KNRM,
			BOLD, filename_pattern, KNRM);
	}
	
	const std::string basename = match[1];
	std::string_view run_n_start_view = {
		name.data() + match.position(3), (size_t)match.length(3)
	};
	Option<u32> maybe_n_start = mnd::stou(run_n_start_view);
	Option<u32> maybe_n_end   = maybe_n_start;
	if(maybe_n_start.is_none()) {
        MND_THROW("mnd::fs:file_info: In file name: %*s , "
			"expected `<start>` run number to be "
			"parsable to uint32_t", name.c_str());
	}

	if(match[3].matched) {
		std::string_view run_n_end_view = {
			name.data() + match.position(3), (size_t)match.length(3)
		};
		maybe_n_end = mnd::stou(run_n_end_view);
        if(maybe_n_end.is_none()) {
			MND_THROW("mnd::fs:file_info: In file name: %*s , "
				"expected `<end>` run number to be "
				"parsable to uint32_t", name.c_str());
		}
	}

    return { std::move(basename), maybe_n_start.unwrap(), maybe_n_end.unwrap()};
}

u32 mnd::fs::file_start_number(std::string_view file) {
    return std::get<1>( file_info(file) );
}

u32 mnd::fs::file_end_number(std::string_view file) {
    return std::get<2>( file_info(file) );
}

std::pair<u32, u32>
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

    return std::string{file_prefix}
        + "_"
        + mnd::utos(bounds.first, nchars_run_number)
        + "_"
        + mnd::utos(bounds.second, nchars_run_number);
}
