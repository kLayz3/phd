#pragma once

#include <optional>
#include <functional>

/* Rust fanboi? Hell yea! Simple:
 * enum Option<T> {
 *   Some(T),
 *   None
 * }
 * A true algebraic sum type! Nullability isn't tied to
 * a self-defined `nil` subset within `T` itself.
 *
 * ~~ Keep wearing your programming socks ~~~
 * ⣇⣿⠘⣿⣿⣿⡿⡿⣟⣟⢟⢟⢝⠵⡝⣿⡿⢂⣼⣿⣷⣌⠩⡫⡻⣝⠹⢿⣿⣷
 * ⡆⣿⣆⠱⣝⡵⣝⢅⠙⣿⢕⢕⢕⢕⢝⣥⢒⠅⣿⣿⣿⡿⣳⣌⠪⡪⣡⢑⢝⣇
 * ⡆⣿⣿⣦⠹⣳⣳⣕⢅⠈⢗⢕⢕⢕⢕⢕⢈⢆⠟⠋⠉⠁⠉⠉⠁⠈⠼⢐⢕⢽
 * ⡗⢰⣶⣶⣦⣝⢝⢕⢕⠅⡆⢕⢕⢕⢕⢕⣴⠏⣠⡶⠛⡉⡉⡛⢶⣦⡀⠐⣕⢕
 * ⡝⡄⢻⢟⣿⣿⣷⣕⣕⣅⣿⣔⣕⣵⣵⣿⣿⢠⣿⢠⣮⡈⣌⠨⠅⠹⣷⡀⢱⢕
 * ⡝⡵⠟⠈⢀⣀⣀⡀⠉⢿⣿⣿⣿⣿⣿⣿⣿⣼⣿⢈⡋⠴⢿⡟⣡⡇⣿⡇⡀⢕
 * ⡝⠁⣠⣾⠟⡉⡉⡉⠻⣦⣻⣿⣿⣿⣿⣿⣿⣿⣿⣧⠸⣿⣦⣥⣿⡇⡿⣰⢗⢄
 * ⠁⢰⣿⡏⣴⣌⠈⣌⠡⠈⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣬⣉⣉⣁⣄⢖⢕⢕⢕
 * ⡀⢻⣿⡇⢙⠁⠴⢿⡟⣡⡆⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣵⣵⣿
 * ⡻⣄⣻⣿⣌⠘⢿⣷⣥⣿⠇⣿⣿⣿⣿⣿⣿⠛⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
 * ⣷⢄⠻⣿⣟⠿⠦⠍⠉⣡⣾⣿⣿⣿⣿⣿⣿⢸⣿⣦⠙⣿⣿⣿⣿⣿⣿⣿⣿⠟
 * ⡕⡑⣑⣈⣻⢗⢟⢞⢝⣻⣿⣿⣿⣿⣿⣿⣿⠸⣿⠿⠃⣿⣿⣿⣿⣿⣿⡿⠁⣠
 * ⡝⡵⡈⢟⢕⢕⢕⢕⣵⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣶⣿⣿⣿⣿⣿⠿⠋⣀⣈⠙
 * ⡝⡵⡕⡀⠑⠳⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠛⢉⡠⡲⡫⡪⡪⡣
 *
 * This type *should not* be used as direct replacement of STL's optional.
 * Is used mainly for monadic operations to have nice representation pre-C++20.
 * Alias for the true std::optional in MONAD is `Maybe` type.
 * Reason is that various perfect forwarding, correct reference returning
 * machinery isn't fully tested or built in. */
/* NOTE: No deref/bool operator, unlike STL optional. This is to keep the intent very explicit in the code. */

namespace mnd {
namespace type_traits {
#if __cplusplus >= 202002L /* Mirrors STL-terminology */

template<typename T>
using remove_cvref = std::remove_cvref<T>;
template<typename T>
using remove_cvref_t = std::remove_cvref_t<T>;

#else

template<typename T>
struct remove_cvref {
	using type = std::remove_cv_t<std::remove_reference_t<T>>;
};
template<typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

#endif // __cplusplus >= 202002l

template<typename L, typename R, typename = void>
struct is_comparable : std::false_type {};

template<typename L, typename R>
struct is_comparable<L, R, std::void_t<
	decltype(std::declval<L>() == std::declval<R>())
>> : std::is_convertible<
	decltype(std::declval<L>() == std::declval<R>()),
	bool
> {};

} // namespace type_traits

using None_t = std::nullopt_t;
inline constexpr auto None = std::nullopt;

template<typename T>
struct Some {
	T value;
};
template<typename T>
Some(T) -> Some<T>;

template<typename T>
class Option {
public:
	using value_type = T; // needed for CLI11

	constexpr Option() : data(None) {};
	constexpr Option(std::nullopt_t) : data(None) {};
	
	template<typename U>
	constexpr Option(Some<U> some) : data( T{ std::move(some.value) } ) {}

	constexpr bool is_some() const noexcept { return data.has_value(); }
	constexpr bool is_none() const noexcept { return !is_some(); }

	/* May panic (throw). Unlike rust, returns back a reference when called on lvalue. */
	constexpr const T&  unwrap() const&  { return data.value(); }
	constexpr T&        unwrap() &       { return data.value(); }
	constexpr T&&       unwrap() &&      { return std::move(data.value()); }
	constexpr const T&& unwrap() const&& { return std::move(data.value()); }

	constexpr decltype(auto) get() const noexcept { return (data); }
	constexpr decltype(auto) get() noexcept { return (data); }

	/* Normally in STL, the functor type `F` is constrained by different concepts. */	
	template<typename U = std::remove_cv_t<T>>
	constexpr T value_or( U&& default_value ) const& {
		return is_some() ? unwrap() : static_cast<T>(std::forward<U>(default_value));
	}
	template<typename U = std::remove_cv_t<T>>
	constexpr T value_or( U&& default_value ) && {
		return is_some() ? std::move(unwrap()) : static_cast<T>(std::forward<U>(default_value));
	}
	template<typename F>
	constexpr auto and_then(F&& f) & {
		return is_some()
			? std::invoke(std::forward<F>(f), data.value())
			: type_traits::remove_cvref_t<std::invoke_result_t<F, T&>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const& {
		return is_some()
			? std::invoke(std::forward<F>(f), data.value())
			: type_traits::remove_cvref_t<std::invoke_result_t<F, T const&>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) && {
		return is_some()
			? std::invoke(std::forward<F>(f), std::move(data.value()))
			: type_traits::remove_cvref_t<std::invoke_result_t<F, T>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const&& {
		return is_some()
			? std::invoke(std::forward<F>(f), std::move(data.value()))
			: type_traits::remove_cvref_t<std::invoke_result_t<F, T const>>{};
	}
	template<typename F>
	constexpr Option or_else( F&& f ) const& {
		return is_some() ? *this : std::forward<F>(f)();
	};
	template<typename F>
	constexpr Option or_else( F&& f ) && {
		return is_some() ? std::move(*this) : std::forward<F>(f)();
	};

	/* fn map<U, F>(self, f: F) -> Option<U> */
	template<typename F>
	constexpr auto map(F&& f) & {
		using U = std::remove_cv_t<std::invoke_result_t<F, T&>>;
		return is_some()
			? Option<U>{Some<U>{std::invoke(std::forward<F>(f), data.value())}}
			: Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) const& {
		using U = std::remove_cv_t<std::invoke_result_t<F, const T&>>;
		return is_some()
			? Option<U>{Some<U>{std::invoke(std::forward<F>(f), data.value())}}
			: Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) && {
		using U = std::remove_cv_t<std::invoke_result_t<F, T&&>>;
		return is_some()
			? Option<U>{Some<U>{std::invoke(std::forward<F>(f), std::move(data).value())}}
			: Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) const&& {
		using U = std::remove_cv_t<std::invoke_result_t<F, const T&&>>;
		return is_some()
			? Option<U>{Some<U>{std::invoke(std::forward<F>(f), std::move(data).value())}}
			: Option<U>{None};
	}

	/* Reset the state back to the `No` variant. */
	constexpr void reset() noexcept { data.reset(); }
	
	template<typename U>
	Option& operator=(U&& rhs) {
		data = std::forward<U>(rhs);
		return *this;
	}
	template<typename U>
	Option& operator=(Some<U>&& rhs) {
		data = std::move( rhs.value );
		return *this;
	}
	template<typename U>
	Option& operator=(const Some<U>& rhs) {
		data = rhs.value;
		return *this;
	}

	/* Explicit comparison only allowed against another Option<U> types :) */
	template<typename U,
		typename = std::enable_if_t<
			type_traits::is_comparable<const T&, const U&>::value
	>> constexpr bool operator==(const Option<U>& rhs) const {
		if(is_some() && rhs.is_some())
			return unwrap() == rhs.unwrap();

		return is_none() == rhs.is_none();
	}
	template<typename U,
		typename = std::enable_if_t<
			type_traits::is_comparable<const T&, const U&>::value
	>> constexpr bool operator!=(const Option<U>& rhs) const {
		return !(*this == rhs);
	}

	constexpr bool operator==(None_t ) const noexcept {
		return is_none();
	}
	constexpr bool operator!=(None_t ) const noexcept {
		return is_some();
	}
	friend constexpr bool operator!=(None_t , const Option& rhs) noexcept {
		return rhs.is_some();
	}
	friend constexpr bool operator==(None_t , const Option& rhs) noexcept {
		return rhs.is_none();
	}
	/* ^^^ all possible comparisons given. */

protected:
	std::optional<T> data;
};

} // namespace mnd
