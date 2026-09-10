#include "FromChars.h"

#include <iterator>
#include <sstream>
#include <charconv>
#include <iomanip>

using i32 =  int32_t;
using u32 = uint32_t;

using namespace mnd;

template <
	typename T,
	bool MatchWholeView = true,
	bool MunchView = false
> static Option<T> string_to_i_helper_(std::string_view& sv) {
	if(sv.empty())
		return None;

	T value{};
	auto [ptr, ec] = std::from_chars(
		sv.data(),
		sv.data() + sv.size(),
		value
	);

	if(ec != std::errc{}) {
		return None;
	}

	if constexpr(MatchWholeView) {
		 if( ptr != sv.data() + sv.size() )
			return None;
	}
	if constexpr(MunchView) {
		sv.remove_prefix(
			std::distance(sv.data(), ptr)
		);
	}

	return Some{value};
}

template<typename T>
static std::string i_to_string_helper_(
	T value,
	u32 width,
	char fill_char
) {
	std::ostringstream oss;
	oss << std::internal
		<< std::setw( static_cast<int>(width) )
		<< std::setfill(fill_char)
		<< value;

	return oss.str();
}

mnd::Option<i32> mnd::stoi(std::string_view sv) {
	return string_to_i_helper_<i32>(sv);
}
mnd::Option<u32> mnd::stou(std::string_view sv) {
	return string_to_i_helper_<u32>(sv);
}

std::string mnd::itos(i32 value, u32 width, char fill_char) {
	return i_to_string_helper_(value, width, fill_char);
}
std::string mnd::utos(u32 value, u32 width, char fill_char) {
	return i_to_string_helper_(value, width, fill_char);
}

mnd::Option<i32> mnd::stoi_munch(std::string_view& sv) {
	return string_to_i_helper_<i32, false, true>(sv);
}
mnd::Option<u32> mnd::stou_munch(std::string_view& sv) {
	return string_to_i_helper_<u32, false, true>(sv);
}
