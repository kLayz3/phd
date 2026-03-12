#pragma once

#include <array>

namespace mnd {

template<typename T, typename U>
bool IsInside(const T& value, const std::array<U, 2>& bounds) {
	static_assert (
        std::is_convertible_v<decltype(std::declval<const U&>() <= std::declval<const T&>()), bool> &&
        std::is_convertible_v<decltype(std::declval<const T&>() <  std::declval<const U&>()), bool>,
        "T and U must support comparison operators T <= U and U < T."
    );
	return bounds[0] <= value and value < bounds[1];
}

}
