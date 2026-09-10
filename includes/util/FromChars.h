#pragma once

#include <cstdint>
#include <string_view>
#include <string>

#include "util/Option.hxx"

namespace mnd {

Option<int32_t> stoi(std::string_view );
Option<uint32_t> stou(std::string_view );

std::string itos(int32_t, uint32_t, char = '0');
std::string utos(uint32_t, uint32_t, char = '0');

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
