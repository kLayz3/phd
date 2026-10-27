#include "FromChars.h"

#include <charconv>
#include <cmath>
#include <limits>
#include <locale>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <type_traits>

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

template <
	typename T,
	bool MatchWholeView = true,
	bool MunchView = false
> static Option<T> string_to_f_helper_(std::string_view& sv_) {
	static_assert(std::is_floating_point_v<T>, "Type must be floating point.");
	if(sv_.empty())
		return None;
	
	std::string_view sv = sv_;
	T value{};
	
	/* Accept +1.0233, e.g. test if plus is upfront and if a digit follows. */
	if(sv[0] == '+') {
		if(sv.size() < 2 || sv[1] == '+' || sv[1] == '-')
			return None;

		sv.remove_prefix(1);
		
		auto parsed = string_to_f_helper_<T, MatchWholeView, MunchView>(sv);
		if constexpr(MunchView) {
			if(parsed.has_value())
				sv_ = sv;
		}

		return parsed;
	}

	const auto [ptr, ec] = std::from_chars(
		sv.data(), sv.data() + sv.size(), value
	);
	if(ec != std::errc{} && ec != std::errc::result_out_of_range)
		return None;

	if constexpr(MatchWholeView) {
		if(ptr != sv.data() + sv.size())
			return None;
	}

	if(ec == std::errc::result_out_of_range) {
		// `from_chars` leaves value unchanged on both overflow and underflow.
		// C-locale distinguishes them w/o changing them and can also recover infs/nans.
		std::istringstream stream{std::string{sv.data(), ptr}};
		stream.imbue(std::locale::classic());
		stream >> value;
		if(stream.fail()) {
			if(value == std::numeric_limits<T>::max()
			   || value == -std::numeric_limits<T>::max()) {
				value = std::copysign(std::numeric_limits<T>::infinity(), value);
			} else {
				return None;
			}
		}
	}

	if constexpr(MunchView) {
		sv.remove_prefix(static_cast<size_t>(ptr - sv.data()));
		sv_ = sv;
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
template<typename T>
static std::string f_to_string_helper_(
	T value,
	u32 precision,
	u32 width,
	char fill_char,
	std::chars_format fmt
) {
	using U = mnd::type_traits::remove_cvref_t<T>;
	static_assert(std::is_floating_point_v<U>, "Type passed must be a floating point.");

	std::ostringstream oss;
	switch(fmt) {
		case std::chars_format::fixed:
			oss << std::fixed;
			break;
		case std::chars_format::scientific:
			oss << std::scientific;
			break;
		case std::chars_format::hex:
			oss << std::hexfloat;
			break;
		case std::chars_format::general:
			oss << std::defaultfloat;
			break;
		default:
			throw std::invalid_argument("Invalid floating-point format passed to mnd::ftos");
	}

	oss << std::internal
		<< std::setprecision( static_cast<int>(precision) )
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
template<typename T>
std::string mnd::ftos(
	T value,
	u32 precision,
	u32 width,
	char fill_char,
	std::chars_format fmt
) {
	return f_to_string_helper_(value, precision, width, fill_char, fmt);
}
template std::string mnd::ftos<long double>(long double, u32, u32, char, std::chars_format);
template std::string mnd::ftos<double>(double, u32, u32, char, std::chars_format);
template std::string mnd::ftos<float>(float, u32, u32, char, std::chars_format);

template<typename T>
mnd::Option<T> mnd::stof(std::string_view sv) {
	return string_to_f_helper_<T>(sv);
}

template Option<long double> mnd::stof<long double>(std::string_view );
template Option<double> mnd::stof<double>(std::string_view );
template Option<float> mnd::stof<float>(std::string_view );

mnd::Option<i32> mnd::stoi_munch(std::string_view& sv) {
	return string_to_i_helper_<i32, false, true>(sv);
}
mnd::Option<u32> mnd::stou_munch(std::string_view& sv) {
	return string_to_i_helper_<u32, false, true>(sv);
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

