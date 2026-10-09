#include "util/MacroHelpers.h"
#include "util/FromChars.h"

#include <exception>
#include <regex.h>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <regex>
#include <string_view>
#include <sys/wait.h>
#include <unistd.h>

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
	if(!pad) return out;

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

std::vector<std::string> mnd::parse_file_to_lines(const std::string& fileName) {
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

std::string mnd::parse_file_to_string(const std::string& fileName) {
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

std::string mnd::longest_common_prefix(mnd::span<const std::string> sequence) {
	if(sequence.empty())
		return {};

	std::string_view prefix = sequence.front();
	for(const auto& s : sequence.subspan(1)) {
		const auto limit = std::min(prefix.size(), s.size());
		size_t i = 0;

		while(i < limit && prefix[i] == s[i])
			++i;

		prefix = prefix.substr(0, i);
		if(prefix.empty())
			break;
	}

	return std::string{prefix};
}

static bool is_word(std::string_view text) {
	static const std::regex re{R"(\w+)"};
	return std::regex_match(text.begin(), text.end(), re);
}

/* Returns the position of where the corresponding closing brace should be.
 * NOTE that the text view must be left at the exact placement of the bracket!
 * If the closing bracket is not found. Returns std::string_view::npos. */
static size_t find_balanced_brace(std::string_view text) {
	if(text.size() < 2)
		MND_THROW("Text is: '%.*s', with not enough size to have two matching brace?",
			(int)text.size(), text.data());

	const char O = text.front();
	if(O != '(' && O != '[' && O != '{' && O != '<')
		MND_THROW("View isn't placed at one of the opening braces ([{<. It is: '%c'", text.front());
	
	const char C = O == '(' ? ')' :
	               O == '[' ? ']' :
	               O == '{' ? '}' : '>';

	size_t depth = 1;
	/* Problem is, just stopping parsing on first '}' would be a mistake,
	 * as args could have a bracket couple '{}' in them. Either way, the brackets must
	 * be balanced. */
	for(size_t i = 1; i < text.size(); ++i) {
		if(text[i] == O) {
			++depth;
		}
		/* Handle quotes concretely.. */
		else if(text[i] == '\'' || text[i] == '"') {
			const char quote = text[i++];
			for(; i < text.size(); ++i) {
				if(text[i] == quote)
					break;
				if(quote == '"' && text[i] == '\\' && i+1 < text.size())
					++i;
			}
			if(i == text.size())
				return std::string_view::npos;
		}
		else if (text[i] == C) {
			--depth;
			if(depth == 0)
				return i;
		}
	}
	return std::string_view::npos;
}

mnd::Option<std::string_view> mnd::extract_text_body (
	std::string_view block_keyword, // SECTION
	std::string_view label, // section name
	std::string_view text   // full text
) {

	if(!is_word(block_keyword))
		MND_THROW("The label is '%.*s', but should only contain [a-zA-Z0-9_] symbols.",
			(int)block_keyword.size(), block_keyword.data());

	if(!label.empty() && !is_word(label))
		MND_THROW("The label is '%.*s', but should only contain [a-zA-Z0-9_] symbols.",
			(int)label.size(), label.data());

	std::regex re;
	try {
		std::string rgx {}; rgx.reserve(block_keyword.size() + label.size() + 16);
		rgx += R"(\b)";
		rgx += block_keyword;
		rgx += R"(\s*\(\s*)";
		rgx += label;
		rgx += R"(\s*\))";

		re = std::regex{rgx};
	} catch(const std::exception& e) {
		MND_THROW("Regex creation failed: reason: %s", e.what());
	}

	mnd::Option<std::string_view> result;
	
	std::match_results<std::string_view::const_iterator> m;
	if(!std::regex_search(text.begin(), text.end(), m, re))
		return None;
	
	const size_t header_end =
		static_cast<size_t>(m[0].second - text.begin());
	const size_t brace = text.find_first_not_of(" \t\r\n\f\v", header_end);

	if(brace == std::string_view::npos || text[brace] != '{')
		MND_THROW("Expected '{' after section header.");

	text.remove_prefix(brace);
	size_t closing_brace = find_balanced_brace(text);

	if(closing_brace == std::string_view::npos)
		MND_THROW("mnd::extract_text_body: "
			"Parse error: unterminated block-section %s(%.*s)\n",
			(int)block_keyword.size(), block_keyword.data(), (int)label.size(), label.data());
	
	result = Some{ text.substr(1, closing_brace-1) };

	text = text.substr(closing_brace+1);

	if(std::regex_search(text.begin(), text.end(), m, re))
		MND_THROW("Duplicate instance of '%s(%.*s)' block found.\n",
			(int)block_keyword.size(), block_keyword.data(), (int)label.size(), label.data());
	
	return result;
}

std::vector<std::string_view> mnd::extract_available_sections (
	std::string_view block_keyword,
	std::string_view text
) {
	if(!is_word(block_keyword))
		MND_THROW("The label is '%.*s', but should only contain [a-zA-Z0-9_] symbols.",
			(int)block_keyword.size(), block_keyword.data());

	std::vector<std::string_view> results{};

	const auto re = std::regex{
		std::string{R"(\b)"}
		+ std::string{block_keyword}
		+ R"(\s*\(\s*(\w*)\s*\)\s*\{)"
	};
	std::match_results<std::string_view::const_iterator> m;
	
	while(std::regex_search(text.begin(), text.end(), m, re)) {
		const auto start = static_cast<size_t>(m[1].first - text.begin());
		results.push_back(text.substr(
			start, static_cast<size_t>(m[1].length())
		));

		const size_t brace = static_cast<size_t>(m[0].second - text.begin()) - 1;

		// The match includes '{'; position the view on it.
		text.remove_prefix(brace);
		const size_t closing_brace = find_balanced_brace(text);

		if(closing_brace == std::string_view::npos)
			MND_THROW("Unterminated section body %.*s.",
				(int)results.back().size(), results.back().data());

		text.remove_prefix(closing_brace + 1);
	}
	return results;
}

/* https://stackoverflow.com/a/7408245/4487530 */
std::vector<std::string> mnd::split(const std::string &text, char sep) {
	std::vector<std::string> tokens;
	std::string::size_type start = 0, end = 0;
	while((end = text.find(sep, start)) != std::string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}

std::vector<std::string_view> mnd::to_views(const std::vector<std::string>& seq) {
	std::vector<std::string_view> views;
	views.reserve(seq.size());
	for(const std::string& str : seq)
		views.emplace_back(str);
	return views;
}

std::vector<std::string_view> mnd::split_view(std::string_view text, char sep) {
	std::vector<std::string_view> tokens;
	std::string::size_type start = 0, end = 0;
	while((end = text.find(sep, start)) != std::string::npos) {
		tokens.emplace_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.emplace_back(text.substr(start));
	return tokens;
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

    if(name.size() > file.size())
		MND_THROW("Invalid file path? %.*s", (int)file.size(), file.data());

	std::smatch match;
    if(!std::regex_match(name, match, re)) {
		MND_THROW("mnd::fs::file_info(): provided file name '%s' does not "
			"match the regular expression: %s%s%s\n",
			name.c_str(), BOLD, filename_pattern, KNRM);
	}
	
	const std::string basename = match[1];
	std::string_view run_n_start_view = {
		name.data() + match.position(2), (size_t)match.length(2)
	};
	const Option<u32> maybe_n_start = mnd::stou(run_n_start_view);
	Option<u32> maybe_n_end   = maybe_n_start;
	if(maybe_n_start.is_none()) {
        MND_THROW("mnd::fs:file_info: In file name: %s , "
			"expected `<start>` run number to be "
			"parsable to uint32_t", name.c_str());
	}

	if(match[3].matched) {
		std::string_view run_n_end_view = {
			name.data() + match.position(3), (size_t)match.length(3)
		};
		maybe_n_end = mnd::stou(run_n_end_view);
        if(maybe_n_end.is_none()) {
			MND_THROW("mnd::fs:file_info: In file name: %s , "
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

std::filesystem::path mnd::fs::resolve_maybe_symlink(const std::filesystem::path& path) {
	if(!std::filesystem::is_symlink(path))
		return path;

	auto target = std::filesystem::read_symlink(path);
	return target.is_absolute() ? target : path.parent_path() / target;
}
