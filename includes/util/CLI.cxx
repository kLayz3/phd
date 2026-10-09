#include "CLI.h"
#include <cstring>
#include <string>

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


#ifdef __linux__

#include <unistd.h>

std::filesystem::path mnd::fs::current_executable_path() {
	std::array<char, 4096> buf{};

	ssize_t n = ::readlink("/proc/self/exe", buf.data(), buf.size() - 1);
	if(n < 0)
		MND_THROW("readlink(/proc/self/exe) failed");

	buf[ (size_t)n ] = '\0';
	return std::filesystem::path(buf.data());
}
std::filesystem::path mnd::fs::current_executable_name() {
	return current_executable_path().stem();	
}

#endif // __linux__
