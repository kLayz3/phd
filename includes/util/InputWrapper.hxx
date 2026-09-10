#pragma once

#include <utility>
#include "CLI.h"

namespace mnd {

/* Small little wrapper to allow strong typing.
 * This is basically a zero-cost abstraction that allows really
 * pretty API's to directly name the positional arguments. */
template<typename T, typename Tag = void>
struct InputWrapper {
	using value_type = T;
	T value;
	InputWrapper() = default;
	InputWrapper(T v) : value(std::move(v)) {}

	operator T&() noexcept { return value; }
	operator const T&() const noexcept { return value; }

    T&       get() &       noexcept { return value; }
    T const& get() const & noexcept { return value; }
    T&&      get() &&      noexcept { return std::move(value); }
};

/* Custom char buffer streaming operations for the phantom wrapper types, if the underlying type
 * implements them. If underlying type's definitions are not found at this point, then this
 * template is sfinae'd out. E.g. vector|array overload is in `json_struct_def.hh`, and won't be
 * automatically detected here, if that specific header is included *after* this one.
 * ADL resolves the namespace mnd:: here.
 *
 * Non-templated specialized overloads can still be defined and compiler will like them more. Obviously. */
template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::type_traits::is_istreamable<T>::value>
> std::istream& operator>>(std::istream& in, InputWrapper<T, Tag>& value) {
    return in >> value.get();
}

template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::type_traits::is_ostreamable<T>::value>
> std::ostream& operator<<(std::ostream& out, InputWrapper<T, Tag> const& value) {
    return out << value.get();
}

} // namespace mnd
