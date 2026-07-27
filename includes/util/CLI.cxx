#include "CLI.h"
#include <cstring>
#include <string>
#include <type_traits>

mnd::Maybe<std::string_view> mnd::extract_text_body (
		std::string_view block_name,
		std::string_view text, 
		std::string_view label
) {
	std::string needle;
	needle.reserve(block_name.size() + label.size() + 2);
	needle += block_name;
	needle += "(";
	needle += label;
	needle += ")";

	std::size_t pos = 0;

	while(true) {
		pos = text.find(needle, pos);
		if(pos == std::string_view::npos)
			return mnd::None;

		std::size_t brace = text.find('{', pos + needle.size());
		if(brace == std::string_view::npos)
			return std::nullopt;

		size_t body_start = brace + 1; // where the slice will start
		size_t depth = 1;

		for(size_t i = brace+1; i < text.size(); ++i) {
			if(text[i] == '{') {
				++depth;
			} else if (text[i] == '}') {
				--depth;

				if(depth == 0) {
					return text.substr(body_start, i - body_start);
				}
			}
		}

		ERROR("unterminated section \'%s(%s)\'\n", std::string(block_name).c_str(), std::string(label).c_str());
	}
}

static auto is_space = [](unsigned char c) {
	return std::isspace(c) != 0;
};

/* Main problem is that `storage` vector's underlying buffer
 * could potentially move, and because of SSO (small string optimization)
 * would leave the argv[i] pointers dangling. Therefore, we absolutely reserve
 * the capacity up front, and not never be able to 'resize' the vector. */
int mnd::Argv::push_back(std::string s) {
	if(storage.size() == capacity)
		MND_THROW("Attempting to add an argument \'%s\' to mnd::Argv but capacity %zu is reached.", s.c_str(), capacity);

	storage.push_back( std::move(s) );
	argv.back() = storage.back().data();
	argv.push_back(nullptr); // last entry must be terminated.

	return static_cast<int>( storage.size() );
}

mnd::Argv mnd::parse_argv(std::string_view text, std::string program_name) {
	static constexpr size_t N_MAX_ARGS = 50;
	Argv out(N_MAX_ARGS);
	out.push_back( std::move(program_name) );

	size_t i = 0;

	while(true) {
		while(i < text.size() && is_space(static_cast<unsigned char>(text[i])))
			++i;

		if(i == text.size())
			break;

		std::string arg;

		while(i < text.size() && !is_space(static_cast<unsigned char>(text[i]))) {
			switch(text[i]) {
				case '\\': { // escape char, just directly dump the raw following character..
					++i;
					if(i == text.size()) ERROR("Dangling backslash in argument list\n");
					arg.push_back(text[i++]);
					break;
				};

				case '\'': { 
					++i;
					while(i < text.size() && text[i] != '\'')
						arg.push_back(text[i++]);

					if(i == text.size()) ERROR("unterminated single quote\n"); 
					// At this point: text[i] == '\'';
					++i;
					break;
				};

				case '"': {
					++i;
					while(i < text.size() && text[i] != '"') {
						if(text[i] == '\\') {
							++i; 
							if(i == text.size()) ERROR("dangling backslash in double quote\n");
							// Will directly encode the i+1 char.
						}

						arg.push_back(text[i++]);
					}

					if(i == text.size())
						ERROR("unterminated double quote\n");

					++i;
				}
				default:
					arg.push_back(text[i++]);
			}
		}
		
		out.push_back( std::move(arg) );
	}
	return out;
	// Now, out.argv[0] and out.argv.data() and is ready to be passed to execve. */
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

bool mnd::parse(std::string_view v, bool& b) {
	if(v == "true" || v == "True" || v == "TRUE" || v == "1") {
		b = true;
		return true;
	}
	if(v == "false" || v == "False" || v == "FALSE" || v == "0") {
		b = false;
		return true;
	}
	return false;
}

#ifdef __linux__

#include <unistd.h>

std::filesystem::path mnd::current_executable_path() {
	std::array<char, 4096> buf{};

	ssize_t n = ::readlink("/proc/self/exe", buf.data(), buf.size() - 1);
	if(n < 0)
		MND_THROW("readlink(/proc/self/exe) failed");

	buf[ (size_t)n ] = '\0';
	return std::filesystem::path(buf.data());
}

#endif
