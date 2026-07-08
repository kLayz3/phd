#pragma once
#include "json_struct_def.hh" // for std::ostream& operator<< on arrays/ranges.
#include "../magic_enum/magic_enum.hpp"
#include "../monad/monad.hxx"

#include <string_view>
#include <charconv>
#include <type_traits>

#define CLI11_ENABLE_EXTRA_VALIDATORS 1
#include "../cli/CLI11.hpp"

using DisplayDefault = mnd::BinaryOpt;

template <
	DisplayDefault d = DisplayDefault::Yes,
	typename T
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	T& variable,
	const std::string& description
) {
	auto* opt = app.add_option(name, variable, description)
		->each ( 
			[name](const std::string& match) {
				WARN("Parsed option "); 
				std::cerr << KBH_YEL << name << KNRM << " as " 
					<< KBH_YEL << match << KNRM << '\n';
			}
		);
	if constexpr(d == DisplayDefault::Yes)	
		opt->capture_default_str();
	return opt;
}

inline CLI::Option* add_logged_flag (
	CLI::App& app,
	std::string name,
	bool& variable,
	std::string description
) {
	return app.add_flag(name, variable, description)
		->each (
			[name](const std::string& match) {
				WARN("Parsed flag "); 
				std::cerr << KBH_YEL << name << KNRM << " as " 
					<< KBH_YEL << match << KNRM << '\n';
			}
		);
}

template <
	DisplayDefault d = DisplayDefault::Yes,
	typename E
> CLI::Option* add_enum_option (
	CLI::App& app,
	const std::string& name,
	E& variable,
	const std::string& description 
) {
	static_assert(std::is_enum_v<E>);

	auto* opt = app.add_option_function<std::string>(
		name, 
		[&variable, &name](const std::string& s) {
			auto e = magic_enum::enum_cast<E>(s);
			if(!e)
				ERROR("Validation error for enum: \'%s\', "
					"passed in \'%s\' which is not parsable.", 
					mnd::type_name<E>().c_str(), s.c_str());
			variable = *e;
		},
		mnd::sstrcat("Enum: ",  magic_enum::enum_names<E>(), ". ", description)
	)
	->each ( 
		[name](const std::string& match) {
			WARN("Parsed option "); 
			std::cerr << KBH_YEL << name << KNRM << " as " 
				<< KBH_YEL << match << KNRM << '\n';
		}
	);
	if constexpr(d == DisplayDefault::Yes)
		opt->default_str(std::string{magic_enum::enum_name(variable)});

	return opt;
}

namespace mnd {

/* A small wrapper to parse out the sections in the config file block. */
Maybe<std::string_view> extract_section_body(std::string_view , std::string_view );

/* Split a string into smaller substrings. */
std::vector<std::string> split(const std::string& , char );

/* Return a vector of views to the underlying sequence of strings. */
std::vector<std::string_view> to_views(const std::vector<std::string>& );

/* Split a string into smaller substrings, and return a view.
 * The reference could dangle! */
std::vector<std::string_view> split_view(std::string_view , char );
std::vector<std::string_view> split_view(std::string&& , char ) = delete;

template<typename T> 
bool parse(std::string_view v, T& out) {
	static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>,
		"Type T must either be integral or floating point.");
	T value{};

	auto [ptr,ec] = std::from_chars(v.data(), v.data() + v.size(), value);

	if(ec != std::errc{} || ptr != v.data() + v.size()) 
		return false;

	out = value;
	return true;
}
bool parse(std::string_view, bool& );

class Argv {
	std::vector<std::string> storage;
	std::vector<char*> argv;
	const size_t capacity;

public:
	Argv() = delete;
	Argv(size_t n) : capacity(n) {
		assert(n > 0 && "Cannot construct argv without at least one entry (program name).");
		storage.reserve(n);
		argv.reserve(n+1); // fills with nullptr: empty value init of ptr = nullptr 
		argv.push_back( nullptr );
		assert(argv.size() == 1);
	}
	int argc() const noexcept {
		return static_cast<int>(argv.size()) - 1;
	}

	char** data() noexcept {
		return argv.data();
	}
	
	int push_back(std::string );
	const char* operator[](size_t i) const noexcept { return argv[i]; }
};

/* Parse a block of text (1st arg) and the program name (2nd argument) into an
 * argv/argc format, which can then be delegated to exec(3) family: execvp */
Argv parse_argv(std::string_view, std::string );

#ifdef __linux__
	/* Get the executable path which invokes this call. For ROOT macros, use
	 * `gInterpreter` API instead. */
	std::filesystem::path current_executable_path(); 
#endif

/* Parse an array range from a text input by a separator 'c' */
template<unsigned char c, typename T, std::size_t N>
std::istream& operator>>(std::istream& in, std::array<T, N>& out) {
	for(size_t i=0; i<N; ++i) {
		if(i != 0) {
			char sep{};
			if(!(in >> sep) || sep != c) {
				in.setstate(std::ios::failbit);
				return in;
			}
		}
		if(!(in >> out[i])) {
			return in;
		}
	}
	return in;	
}
/* ^^^^^^ We don't put this in global namespace as it would wreak havoc on the 
 * ADL and overload resolver. */
} // namespace mnd

