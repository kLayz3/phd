#pragma once

#include <cstdint>
#include <string_view>
#include <string>

#include "util/Option.hxx"

namespace mnd {

/* Nicer repr of std::stoi, that will parse a character sequence
 * to int/unsigned int (stou), it will match the whole sequence to be parsable.
 * Is currenly accepting only standard C locale, as shown in example [3]
 * [1]: stoi("1234") -> Some{1234}
 * [2]: stoi(" 123") -> None
 * [3]: stoi("123'456") -> None
 * [4]: stoi("123text") -> None
 * [5]: stoi("-12")  -> Some{-12}
 * [6]: stoi("3294967295") -> None
 * [7]: stou("3294967295") -> Some{3294967295}
 * [8]: stou("-12") -> None
 */
Option<int32_t> stoi(std::string_view );
Option<uint32_t> stou(std::string_view );

/* Convert int to string:
 * arg #1: value
 * arg #2: padding spaces, will pad to at least this size
 * arg #3: padding char */
std::string itos(int32_t, uint32_t = 0, char = '0');

/* Convert unsigned int to string:
 * arg #1: value
 * arg #2: padding spaces, will pad to at least this size
 * arg #3: padding char */
std::string utos(uint32_t, uint32_t = 0, char = '0');

/* Case where we wish to consume as many chars as possible.
 * On success, removes the parsed prefix and returns its value.
 * On invalid input or overflow, returns None and leaves the view unchanged. */
Option<int32_t> stoi_munch(std::string_view& );
Option<uint32_t> stou_munch(std::string_view& );

/* Remove leading and trailing whitespaces from a view. Doesn't mutate the underlying buffer,
 * only returns the trimmed view object. */
inline constexpr std::string_view trim(std::string_view sv) noexcept {
	constexpr auto whitespace = " \t\n\r\f\v";

	const auto first = sv.find_first_not_of(whitespace);
	if(first == std::string_view::npos)
		return {};

	const auto last = sv.find_last_not_of(whitespace);
	return sv.substr(first, last - first + 1);
}

}
