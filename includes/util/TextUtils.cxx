#include "util/TextUtils.h"
#include "util/FromChars.h"

#include <regex>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>

#include "monad/accessories/typedefs.hxx"

/* This is pure hatred coding. I absolutely despise the std::string|view API in C++,
 * plus sprinkle some const char* C-api on top. Why not. Junge Alter es sieht ziemlich schlimm aus. */

using namespace std::literals;

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

bool mnd::has_prefix(std::string_view text, std::string_view prefix) {
	return !prefix.empty()
		&& text.size() >= prefix.size()
		&& text.substr(0, prefix.size()) == prefix;
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

	// Opening delim.
	const char O = text.front();
	if(O != '(' && O != '[' && O != '{' && O != '<')
		MND_THROW("View currently isn't placed at one of the opening braces: ([{< . It is placed on: '%c'", text.front());

	// Closing delim
	const char C = O == '(' ? ')' :
	               O == '[' ? ']' :
	               O == '{' ? '}' : '>';

	size_t depth = 1;

	/* Problem is, just stopping parsing on first '}' would be a mistake,
	 * as args could have a bracket couple '{}' in them.
	 * Or they could be escaped '\}'.
	 * Or they could be quoted inside '' or "".
	 * As such, the brackets must be balanced. */
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
		if(text[i] == '\\') { // escape char
			if(++i == text.size())
				MND_THROW("find_balanced_brace: dangling \\ escape char at the edge?");
		}
	}
	return std::string_view::npos;
}

std::vector<std::string> mnd::split(std::string_view text, char sep) {
	const auto views = split_view(text, sep);
	std::vector<std::string> tokens;

	tokens.reserve(views.size());

	/* Note, cannot directly materialize a copy since the view requires some
	 * delimiter parsing out. */
	for(auto view : views) {
		view = mnd::trim(view);
		std::string token; token.reserve(view.size());

		char quote = 0;
		for(size_t i = 0; i < view.size(); ++i) {
			const char c = view[i];
			
			/* Priority takes escape char. */
			if(view[i] == '\\') {
				if(++i == view.size()) MND_THROW("split: dangling escape character.");

				token.push_back( view[i] );
			} else if(quote) {
				if(c == quote)
					quote = 0;
				else
					token.push_back(c);
			} else if(c == '\'' || c == '"') {
				quote = c;
			} else {
				token.push_back(c);
			}
		}
		tokens.push_back(std::move(token));
	}

	return tokens;
}

/* Based on the separator, usually ',' , find the substr inbetween two, keeping
 * in mind that possibly stuff gets escaped, and chars in quote '' and "" cannot constitute
 * a valid separator char. */
std::vector<std::string_view> mnd::split_view(std::string_view text, char sep) {
	std::vector<std::string_view> tokens;
	size_t start = 0;
	char quote = 0;

	for(size_t i = 0; i < text.size(); ++i) {
		const char c = text[i];

		/* Escape char \ escapes everything. Unconditionally. */
		if(c == '\\') {
			if(++i == text.size())
				MND_THROW("split_view: dangling escape character.");
		} else if(quote) {
			if(c == quote) quote = 0;
		} else if(c == '\'' || c == '"') {
			quote = c;
		} else if(c == sep) {
			tokens.push_back(text.substr(start, i - start));
			start = i + 1;
		}
	}
	if(quote)
		MND_THROW("split_view: unterminated quote %c somewhere.", quote);

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

static const std::regex re {mnd::rootfile_pattern};

std::tuple<std::string, u32, u32>
mnd::rootfile_info(std::string_view file) {
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
			name.c_str(), BOLD, rootfile_pattern, KNRM);
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

std::string mnd::run_name(std::string_view file) {
	return std::get<0>( rootfile_info(file) );
}
u32 mnd::rootfile_start_number(std::string_view file) {
    return std::get<1>( rootfile_info(file) );
}
u32 mnd::rootfile_end_number(std::string_view file) {
    return std::get<2>( rootfile_info(file) );
}

std::pair<u32, u32>
mnd::rootfile_number_bounds(const std::vector<std::string>& files) {
    if(files.empty())
        throw std::invalid_argument("Cannot determine bounds of an empty file sequence");

    return {
        rootfile_start_number(files.front()),
        rootfile_end_number(files.back())
    };
}

std::string mnd::rootfile_names_concatenated(const std::vector<std::string>& files) {
    auto bounds = rootfile_number_bounds(files);

    return std::string{file_prefix}
        + "_"
        + mnd::utos(bounds.first, nchars_run_number)
        + "_"
        + mnd::utos(bounds.second, nchars_run_number);
}

std::string mnd::file_stem(std::string_view path) {
	return std::filesystem::path{path}.stem().string();
}

mnd::Option<std::string_view> mnd::sdsl::extract_text_body (
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

std::vector<std::string_view> mnd::sdsl::extract_available_sections (
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

		// The match includes '{'; pin the view on top of it.
		text.remove_prefix(brace);
		const size_t closing_brace = find_balanced_brace(text);

		if(closing_brace == std::string_view::npos)
			MND_THROW("Unterminated section body %.*s.",
				(int)results.back().size(), results.back().data());

		text.remove_prefix(closing_brace + 1);
	}
	return results;
}


