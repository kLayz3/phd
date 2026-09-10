#pragma once
#include "json_struct_def.hh" // for std::ostream& operator<< on ranges.
#include "../magic_enum/magic_enum.hpp"
#include "../monad/monad.hxx"

#include <string_view>
#include <charconv>
#include <type_traits>

#define CLI11_ENABLE_EXTRA_VALIDATORS 1
#include "../cli/CLI11.hpp"

#include "Option.hxx"

using DisplayDefault = mnd::BinaryOpt;

namespace mnd::cli::detail {

struct State {
	bool authoritative_seen = false;
	bool current_is_authoritative = false;
};

inline constexpr char auth_sym = '!';
} //namespace mnd::cli::detail

namespace mnd::type_traits {

template<typename T, typename S = std::istream>
using is_istreamable = CLI::detail::is_istreamable<T, S>;
template<typename T, typename S = std::ostream>
using is_ostreamable = CLI::detail::is_ostreamable<T, S>;

} // namespace mnd::type_traits

/* Overload for non-enum types. */
template <
	DisplayDefault d = DisplayDefault::Yes,
	typename T,
	typename std::enable_if<!std::is_enum_v<T>>::type* = nullptr
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	T& variable,
	const std::string& description
) {
	auto state = std::make_shared<mnd::cli::detail::State>();
	const std::string default_value = CLI::detail::to_string(variable);
	auto* opt = app.add_option_function<T>(
		name, 
		[&variable, name, state](const T& match) {
			if(state->current_is_authoritative) {
				// The ! occurrence overrides everything.
				variable = match;
				state->authoritative_seen = true;
				WARN("Parsed %sauthoritative%s option ", BOLD, KNRM); 
				std::cerr << KBH_YEL << name << KNRM << " as " 
					<< KBH_CYN << match << KNRM << '\n';
			} else if(!state->authoritative_seen) {
				variable = match;
				WARN("Parsed option "); 
				std::cerr << KBH_YEL << name << KNRM << " as " 
					<< KBH_CYN << match << KNRM << '\n';
			}
			state->current_is_authoritative = false;
		}, description)
		->transform( [state](std::string input) -> std::string {
			if(!input.empty() && input.front() == mnd::cli::detail::auth_sym) {
				state->current_is_authoritative = true;
				input.erase(input.begin());
			}
			return input; // RVO
		})
		->trigger_on_parse();

	if constexpr(d == DisplayDefault::Yes)	
		opt->default_str(default_value);
	return opt;
}

template <
	DisplayDefault d = DisplayDefault::Yes,
	typename E,
	typename std::enable_if<std::is_enum_v<E>>::type* = nullptr
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	E& variable,
	const std::string& description
) {
	auto state = std::make_shared<mnd::cli::detail::State>();
	const std::string default_value = std::string{magic_enum::enum_name(variable)};
	auto* opt = app.add_option_function<std::string>(
		name, 
		[&variable, name, state](const std::string& s) {
			auto e = magic_enum::enum_cast<E>(s);
			if(!e)
				ERROR("Validation error for enum: \'%s\', "
					"passed in \'%s\' which is not parsable.",
					mnd::type_name<E>().c_str(), s.c_str());
			if(state->current_is_authoritative) {
				// The ! occurrence overrides everything.
				variable = *e;
				state->authoritative_seen = true;
				WARN("Parsed %sauthoritative%s option ", BOLD, KNRM);
				std::cerr << KBH_YEL << name << KNRM << " as "
					<< KBH_CYN << s << KNRM << '\n';
			} else if(!state->authoritative_seen) {
				variable = *e;
				WARN("Parsed option ");
				std::cerr << KBH_YEL << name << KNRM << " as "
					<< KBH_CYN << s << KNRM << '\n';
			}
			state->current_is_authoritative = false;
		},
		mnd::sstrcat("Enum: ",  magic_enum::enum_names<E>(), ". ", description))
		->transform( [state](std::string input) -> std::string {
			if(!input.empty() && input.front() == mnd::cli::detail::auth_sym) {
				state->current_is_authoritative = true;
				input.erase(input.begin());
			}
			return input; // RVO
		})
		->trigger_on_parse();

	if constexpr(d == DisplayDefault::Yes)
		opt->default_str(default_value);

	return opt;
}

/* A more specialized overload to parse a sequence of enum tokens. */
template <
	DisplayDefault d = DisplayDefault::Yes,
	typename E,
	typename std::enable_if<std::is_enum_v<E>>::type* = nullptr
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	std::vector<E>& variable,
	const std::string& description
) {
	auto state = std::make_shared<mnd::cli::detail::State>();

	std::vector<std::string> default_values;
	default_values.reserve(variable.size());

	for(const auto e : variable)
		default_values.emplace_back(magic_enum::enum_name(e));

	auto* opt = app.add_option_function<std::vector<std::string>>(
		name,
		[&variable, name, state](const std::vector<std::string>& matches) {
			std::vector<E> parsed;
			parsed.reserve(matches.size());

			for(const auto& s : matches) {
				auto e = magic_enum::enum_cast<E>(s);
				if(!e)
					ERROR("Validation error for enum: \'%s\', "
						"passed in \'%s\' which is not parsable.",
						mnd::type_name<E>().c_str(), s.c_str());

				parsed.push_back(*e);
			}

			if(state->current_is_authoritative) {
				// The ! occurrence overrides everything.
				variable = std::move(parsed);
				state->authoritative_seen = true;

				WARN("Parsed %sauthoritative%s option ", BOLD, KNRM);
				std::cerr << KBH_YEL << name << KNRM << " as "
					<< KBH_CYN << matches << KNRM << '\n';
			} else if(!state->authoritative_seen) {
				variable = std::move(parsed);
				WARN("Parsed option ");
				std::cerr << KBH_YEL << name << KNRM << " as "
					<< KBH_CYN << matches << KNRM << '\n';
			}
			state->current_is_authoritative = false;
		},
		mnd::sstrcat("Enum: ", magic_enum::enum_names<E>(), ". ", description))
		->transform([state](std::string input) -> std::string {
			if(!input.empty() &&
			input.front() == mnd::cli::detail::auth_sym)
			{
				state->current_is_authoritative = true;
				input.erase(input.begin());
			}

			return input;
		})
		->trigger_on_parse();

	if constexpr(d == DisplayDefault::Yes)
		opt->default_str(CLI::detail::to_string(default_values));

	return opt;
}

/* For the Option<T> wrapper, also expose a CLI tool template specialization
 * to parse it properly, otherwise boilerplate reeks through the code. */
namespace mnd::cli::detail {

inline constexpr char empty_sym = '@';
inline constexpr char reset_sym = '~';

} // namespace mnd::cli::detail

template <
	DisplayDefault d = DisplayDefault::Yes,
	typename T
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	mnd::Option<T>& variable,
	const std::string& description
) {
	static_assert(
		std::is_default_constructible_v<T>,
		"mnd::Option<T> CLI parsing requires default-constructible T"
	);
	/* God I love undocumented API. So basically the tokenising begins *before*
	 * transformers take place. E.g. `--flag=!@` cannot be parsed for Option<array> ...
	 * Just had to dissect the library like usual. HINT: people please write your docs. */
	 
	/* Inside CLI11.hpp:
	 *   using results_t = std::vector<std::string>;
	 *   using callback_t = std::function<bool(const results_t &)>;
	 */

	auto state = std::make_shared<mnd::cli::detail::State>();
	
	CLI::callback_t callback =
	[&variable, name, state](const CLI::results_t& raw) -> bool {
		if(raw.empty())
			return false;

		/* Copy because we're going to strip the authoritative prefix
		 * from the first vector's element before delegating to CLI11. */
		auto input = raw;
		auto& first = input.front();
	
		/* Problem is that for containers, passing a single flag e.g. '~' will result 
		 * in its results vector being padded by empty strings... */
		auto is_magic = [&input, &first](char symbol) {
			return first.size() == 1 &&
			first.front() == symbol &&
			std::all_of(
				std::next(input.begin()),
				input.end(),
				[](const std::string& s) { return s.empty(); }
			);
		};

		const bool authoritative =
			!first.empty() &&
			first.front() == mnd::cli::detail::auth_sym;

		if(authoritative)
			first.erase(first.begin());

		/* A lone '!' isn't a valid value. */
		if(first.empty())
			return false;

		const bool reset = is_magic(mnd::cli::detail::reset_sym);
		const bool empty = is_magic(mnd::cli::detail::empty_sym);

		T match{}; // local placeholder for value

		/* '@' and '~' both yield T{},
		 * otherwise use CLI11's normal conversion for T. */
		if(!reset && !empty) {
			if(!CLI::detail::lexical_conversion<T,T>(input, match))
				return false;
			/* On successful conversion, keep going. */
		}

		/* '!x' always wins
		 * ordinary x is ignored after any !x
		 * a later !x can replace an earlier !x */
		if(!authoritative && state->authoritative_seen)
			return true;

		if(authoritative)
			state->authoritative_seen = true;

		variable = std::move(match);

		if(reset) {
			variable.reset();

			WARN("Parsed %ssum-type reset%s option ",
				authoritative ?  BOLD "authoritative ": "",
				authoritative ? KNRM : ""
			);

			std::cerr << KBH_YEL << name << KNRM << " as "
				<< MND_RGB_COL(204,102,0) << "none" << KNRM << '\n';

			return true;
		}

		if(authoritative) {
			WARN("Parsed %sauthoritative%s sum-type option ", BOLD, KNRM);
		} else {
			WARN("Parsed sum-type option ");
		}

		std::cerr << KBH_YEL << name << KNRM << " as "
			<< (empty ? MND_RGB_COL(255,102,255) "defaulted value: " : "")
			<< KBH_CYN << variable.unwrap() << KNRM << '\n';

		return true;
	};
	
	auto* opt = app.add_option(
		name,
		std::move(callback),
		description
		 + mnd::msg("\n%s\'%c\'%s flag requests an explicitly empty object (but in the value-given variant), "
			        "\n%s\'%c\'%s requests a reset back to the no-value-given variant.",
		            MND_RGB_COL(255,102,255), mnd::cli::detail::empty_sym, KNRM,
		            MND_RGB_COL(204,102,0  ), mnd::cli::detail::reset_sym, KNRM)
	)
	->type_name(CLI::detail::type_name<T>())
	->type_size(
		1, CLI::detail::type_count<T>::value
	)
	->expected(
		CLI::detail::expected_count<T>::value
	)
	->trigger_on_parse();

	if constexpr(d == DisplayDefault::Yes) {
		opt->default_str(
			variable.map([](const auto& value) {
				return CLI::detail::to_string(value);
			}).value_or("none")
		);
	}

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
					<< KBH_CYN << match << KNRM << '\n';
			}
		);
}

namespace mnd {

/* A small wrapper to parse out the sections in the config file block. */
mnd::Option<std::string_view> extract_text_body(std::string_view , std::string_view , std::string_view );

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
		argv.reserve(n+1); // fills with nullptr
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
namespace fs {
std::filesystem::path current_executable_path();
std::filesystem::path current_executable_name();
}
#endif

/* Parse an array range from a text input by a separator 'c' */
template<unsigned char c, typename Cont, 
	typename std::enable_if<mnd::is_an_array_v<Cont>>::type* = nullptr
> std::istream& operator>>(std::istream& in, Cont& out) {
	constexpr size_t N = mnd::is_an_array<Cont>::size;

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

/* Other overload is specifically for dynamically sized objects */
template<unsigned char c, typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& out) {
	T value;
	
	/* Try fetching an entry. Can immediately fail and be empty.
	 * In this case, just promptly return. */
	if(!(in >> value)) 
		return in;
	
	out.push_back(value);

	while(true) {
		const auto next = in.peek();

		if(next == std::char_traits<char>::eof() ||
		   next != c)
		{
			return in;
		}
		in.get(); // consume it

		if(!(in >> value)) {
			in.setstate(std::ios::failbit); // if the separator is consumed, next token must be a valid value.
			return in;
		}
		out.push_back(std::move(value));
	}
	return in;	
}

} // namespace mnd


namespace CLI {

/* Extra validator wrapper */
template<typename T>
Validator RangeOrEmpty(T min, T max) {
	return Validator {
		[range = Range(min, max)](std::string& input) mutable -> std::string {
			return input == "{}" ? std::string{} : range(input);
		},
		std::string(detail::type_name<T>()) +
			" in [" + std::to_string(min) +
			" - " + std::to_string(max) +
			"] or none"
		,
		"RangeOrEmpty"
	};
}

} // namespace CLI
