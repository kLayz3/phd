#pragma once

#include <cstdint>
#include <string_view>
#include <string>
#include <charconv>

#include "util/Option.hxx"

namespace mnd {

/* Nicer repr of std::stoi, that will parse a character sequence
 * to int/unsigned int (stou), it will match the whole sequence to be parsable.
 * Is currently accepting only standard C locale, as shown in example [3]
 * [1]: stoi("1234") -> Some{1234}
 * [2]: stoi(" 123") -> None
 * [3]: stoi("123'456") -> None
 * [4]: stoi("123text") -> None
 * [5]: stoi("-12")  -> Some{-12}
 * [6]: stoi("3294967295") -> None
 * [7]: stou("3294967295") -> Some{3294967295}
 * [8]: stou("-12") -> None
 */
Option<i32> stoi(std::string_view );
Option<u32> stou(std::string_view );

/* Convert int to string:
 * arg #1: value
 * arg #2: padding spaces, will pad to at least this size
 * arg #3: padding char */
std::string itos(i32, u32 = 0, char = '0');

/* Convert unsigned int to string:
 * arg #1: value
 * arg #2: padding spaces, will pad to at least this size
 * arg #3: padding char */
std::string utos(u32, u32 = 0, char = '0');

/* Convert floating point to string:
 * arg #1: value
 * arg #2: decimal places, will pad to at least this size (from the right of decimal pt)
 * arg #3: at least the total width of the string
 * arg #3: padding char
 * arg #4: formatter type: scientific|fixed|hex|general. Default: general */
template<typename T>
std::string ftos(
	T value,
	u32 precision = 2,
	u32 width = 0,
	char fill_char = '0',
	std::chars_format fmt = std::chars_format::fixed
);

/* Convert string to one of the floats: float|double|long double
 * Is currently accepting only standard C locale, as shown in example [3]
 * E.g., say we deal with `double` type:
 * [ 1]: stof("1234")       -> Some{1234.0}
 * [ 2]: stof("1234.0")     -> Some{1234.0}
 * [ 3]: stof(" 123")       -> None
 * [ 4]: stof("123'456.12") -> None
 * [ 5]: stof("123text")    -> None
 * [ 6]: stof("-12")        -> Some{-12.0}
 * [ 7]: stof("-12.57348")  -> Some{-12.57348}
 * [ 8]: stof("1.02e3")     -> Some{1020.0}
 * [ 9]: stof("1.02 e3")    -> None , because of whitespace
 * [10]: stof("-1.12345E3") -> Some{-1123.45}
 * [11]: stof("1.02e400")   -> Some{std::numeric_limits<double>::infinity} , because of pos. overflow
 * [12]: stof("-1.02e400")  -> Some{-std::numeric_limits<double>::infinity}, because of neg. overflow
 * [13]: stof("NAN")        -> Some{NAN}
 * [14]: stof("inf")        -> Some{std::numeric_limits<double>::infinity}
 * [15]: stof("-inf")       -> Some{-std::numeric_limits<double>::infinity}
 * [16]: stof("+1.02e3")    -> Some{1020.0}
 * [17]: stof("++1.02e3")   -> None , same as other 2-sign combos --,+-,-+ 
 * Is case insensitive, so NAN,nan,inf,INF etc are all accepted. Same for scientific notation e.
 */
template<typename T>
Option<T> stof(std::string_view );

/* Case where we wish to consume as many chars as possible.
 * On success, removes the parsed prefix and returns its value.
 * On invalid input or overflow, returns None and leaves the view unchanged. */
Option<i32> stoi_munch(std::string_view& );
Option<u32> stou_munch(std::string_view& );

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

} // namespace mnd
