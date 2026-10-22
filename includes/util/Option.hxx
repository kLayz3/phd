#pragma once

#include <optional>
#include <utility>
#include <functional>
#include <type_traits>

#include <util/FloatAsserts.hxx>

/* Rust fanboi? Hell yea! Simple:
 * enum Option<T> {
 *   Some(T),
 *   None
 * }
 * A true algebraic sum type! Nullability isn't tied to
 * a self-defined `nil` subset within `T` itself.
 *
 * ~~ Keep wearing your programming socks ~~
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
 * This type *should not* be always used as direct replacement of STL's optional.
 * Is used mainly for monadic operations to have nice representation pre-C++20.
 * Alias for the true std::optional in MONAD is `Maybe` type.
 * Reason is that various perfect forwarding, correct reference returning
 * machinery isn't fully tested or built in. The type sitting nested in the Option<T>
 * type should always be either a non-ref, or a raw pointer. Also non C-array, because obvious.
 *
 * Bigger reason is to have Rust-like niche discriminant packing.
 * Generally in C++ (as of 2026): sizeof(std::optional<T>) > sizeof(T) , because of alignment.
 * However for certain objects such as pointers or non-zero integers, the discriminant can be packed
 * into one of the bit-patterns, thus reducing it to exactly sizeof(T).
 * This is kinda cool to keep the cache playground tidy and packed. */
/* NOTE: No bool operator, unlike STL optional. This is to keep the intent very explicit in the code.
 * NOTE: No direct assignment supported, a'la Rust. `mnd::Option<i32> i = 32;` is a compiler error!
 * Right side needs to be explicitly wrapped as `Some{i32}` .
 * It compiles down efficiently to a single optimised store though. */

namespace mnd {

template<typename T> class Option;

namespace type_traits {
#if __cplusplus >= 202002L /* Mirror STL-terminology */

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

#endif // __cplusplus >= 202002L

template<typename L, typename R, typename = void>
struct is_comparable : std::false_type {};
template<typename L, typename R>
struct is_comparable<L, R, std::void_t<
	decltype(std::declval<L>() == std::declval<R>())
>> : std::is_convertible<
	decltype(std::declval<L>() == std::declval<R>()),
	bool
> {};

template<typename>
struct is_option : std::false_type {};
template<typename T>
struct is_option<Option<T>> : std::true_type {};

} // namespace type_traits

using None_t = std::nullopt_t;
inline constexpr auto None = std::nullopt;

template<typename T>
struct Some {
	T value;
};
template<typename T>
Some(T) -> Some<T>;

/* Abusing CRTP. Lol. Check M. Bancilla, Template Metaprogramming with C++,
 * best invested €38 ever. */
namespace detail {

template<
	typename Derived,
	typename Repr
> class OptionBase {
public:
	using value_type = typename Repr::value_type; // needed for CLI11
	using T = value_type;

	constexpr OptionBase() : data(None) {};
	constexpr OptionBase(None_t) : data(None) {};
	
	template<typename U>
	constexpr OptionBase(Some<U> some) : data( T{ std::move(some.value) } ) {}

	constexpr bool is_some() const noexcept { return data.has_value(); }
	constexpr bool is_none() const noexcept { return !is_some(); }

	/* May panic (throw). Unlike Rust, returns back a reference when called on lvalue. */
	constexpr const T&  unwrap() const&  { return data.value(); }
	constexpr T&        unwrap() &       { return data.value(); }
	constexpr T&&       unwrap() &&      { return std::move(data.value()); }
	constexpr const T&& unwrap() const&& { return std::move(data.value()); }
	
	/* Unchecked access. Delegates to inner's deref op, returns (hopefully) T& category. */
	constexpr decltype(auto) operator*() &       noexcept { return *data; }
	constexpr decltype(auto) operator*() const&  noexcept { return *data; }
	constexpr decltype(auto) operator*() &&      noexcept { return *std::move(data); }
	constexpr decltype(auto) operator*() const&& noexcept { return *std::move(data); }

	/* Normally in STL, the functor type `F` is constrained by different concepts.
	 * This call always returns by value. */
	template<typename U = std::remove_cv_t<T>>
	constexpr T value_or( U&& default_value ) const& {
		return is_some() ? *data : static_cast<T>(std::forward<U>(default_value));
	}
	template<typename U = std::remove_cv_t<T>>
	constexpr T value_or( U&& default_value ) && {
		return is_some() ? *std::move(data) : static_cast<T>(std::forward<U>(default_value));
	}

	/* fn and_then<U, F>(self, f: F) -> Option<U>, where F(T) -> Option<U> */
	template<typename F>
	constexpr auto and_then(F&& f) & {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T&>
		>;
		static_assert(type_traits::is_option<R>::value,
			"mnd::Option<T>::and_then : callback function must return an mnd::Option<U> type"
		);
		return is_some() ? std::invoke(std::forward<F>(f), *data) : R{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const& {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T const&>
		>;
		static_assert(type_traits::is_option<R>::value,
			"mnd::Option<T>::and_then : callback function must return an mnd::Option<U> type"
		);
		return is_some() ? std::invoke(std::forward<F>(f), *data) : R{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) && {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T&&>
		>;
		static_assert(type_traits::is_option<R>::value,
			"mnd::Option<T>::and_then : callback function must return an mnd::Option<U> type"
		);
		return is_some() ? std::invoke(std::forward<F>(f), *std::move(data)) : R{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const&& {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T const&&>
		>;
		static_assert(type_traits::is_option<R>::value,
			"mnd::Option<T>::and_then : callback function must return an mnd::Option<U> type"
		);
		return is_some() ? std::invoke(std::forward<F>(f), *std::move(data)) : R{};
	}

	/* fn or_else<F>(self, f: F) -> Option<T>, where F() -> Option<T> */
	template<typename F>
	constexpr Derived or_else( F&& f ) const& {
		static_assert(std::is_same_v<
			type_traits::remove_cvref_t<std::invoke_result_t<F>>, Derived
		>, "Callback F must return mnd::Option<T>");
		return is_some() ? self() : std::invoke(std::forward<F>(f));
	};
	template<typename F>
	constexpr Derived or_else( F&& f ) && {
		static_assert(std::is_same_v<
			type_traits::remove_cvref_t<std::invoke_result_t<F>>, Derived
		>, "Callback F must return mnd::Option<T>");
		return is_some() ? std::move(self()) : std::invoke(std::forward<F>(f));
	};

	/* fn map<U, F>(self, f: F) -> Option<U>, where F(T) -> U */
	template<typename F>
	constexpr auto map(F&& f) & {
		using U = type_traits::remove_cvref_t<std::invoke_result_t<F, T&>>;
		return is_some()
			? ::mnd::Option<U>{Some<U>{std::invoke(std::forward<F>(f), *data)}}
			: ::mnd::Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) const& {
		using U = type_traits::remove_cvref_t<std::invoke_result_t<F, const T&>>;
		return is_some()
			? ::mnd::Option<U>{Some<U>{std::invoke(std::forward<F>(f), *data)}}
			: ::mnd::Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) && {
		using U = type_traits::remove_cvref_t<std::invoke_result_t<F, T&&>>;
		return is_some()
			? ::mnd::Option<U>{Some<U>{std::invoke(std::forward<F>(f), *std::move(data)) }}
			: ::mnd::Option<U>{None};
	}
	template<typename F>
	constexpr auto map(F&& f) const&& {
		using U = type_traits::remove_cvref_t<std::invoke_result_t<F, const T&&>>;
		return is_some()
			? ::mnd::Option<U>{Some<U>{std::invoke(std::forward<F>(f), *std::move(data)) }}
			: ::mnd::Option<U>{None};
	}

	/* fn filter<F>(self, f: F) -> Option<T>, where F(const T&) -> bool
	 * Returns exact same the same object (by value),
	 * if it passes the predicate, else it returns None. */
	template<typename F>
	constexpr Derived filter(F&& f) const& {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T const&>
		>;
		static_assert(std::is_convertible_v<R, bool>, "mnd::Option<T>::filter : "
			"callback function must accept const T&, and return a type convertible to bool.");

		if(is_some() && std::invoke(std::forward<F>(f), *std::as_const(data)))
			return self();

		return None;
	}
	template<typename F>
	constexpr Derived filter(F&& f) && {
		using R = type_traits::remove_cvref_t<
			std::invoke_result_t<F, T const&>
		>;
		static_assert(std::is_convertible_v<R, bool>, "mnd::Option<T>::filter : "
			"callback function must accept const T&, and return a type convertible to bool.");
			
		if(is_some() && std::invoke(std::forward<F>(f), *std::as_const(data)))
			return std::move(self());

		return None;
	}

	constexpr void reset() noexcept { data.reset(); }
	
	template<typename U>
	constexpr Derived& operator=(Some<U>&& rhs) {
		static_assert(std::is_assignable_v<Repr&, U&&>,
			"Right hand-side inner value type must be move-assignable to the optional representation.");
		data = std::move( rhs.value );
		return self();
	}
	template<typename U>
	constexpr Derived& operator=(const Some<U>& rhs) {
		static_assert(std::is_assignable_v<Repr&, const U&>,
			"Right hand-side inner value type must be copy-assignable to the optional representation.");
		data = rhs.value;
		return self();
	}
	constexpr Derived& operator=(None_t ) noexcept {
		reset();
		return self();
	}

	/* Explicit comparison only allowed against another Option<U> types :) */
	template <
		typename D,
		typename R,
		typename U = typename R::value_type,
		typename = std::enable_if_t<
			type_traits::is_comparable<const T&, const U&>::value
		>
	> constexpr bool operator==(const OptionBase<D,R>& rhs) const {
		if(is_some() && rhs.is_some())
			return **this == *rhs;

		return is_none() == rhs.is_none();
	}
	template <
		typename D,
		typename R,
		typename U = typename R::value_type,
		typename = std::enable_if_t<
			type_traits::is_comparable<const T&, const U&>::value
	>> constexpr bool operator!=(const OptionBase<D,R>& rhs) const {
		return !(*this == rhs);
	}

	constexpr bool operator==(None_t ) const noexcept {
		return is_none();
	}
	constexpr bool operator!=(None_t ) const noexcept {
		return is_some();
	}
	friend constexpr bool operator!=(None_t , const OptionBase& rhs) noexcept {
		return rhs.is_some();
	}
	friend constexpr bool operator==(None_t , const OptionBase& rhs) noexcept {
		return rhs.is_none();
	}

protected:
	Repr data;

private:
	constexpr Derived& self() noexcept {
		return static_cast<Derived&>(*this);
	}
	constexpr const Derived& self() const noexcept {
		return static_cast<const Derived&>(*this);
	}
};

namespace unsafe {

#if defined(__clang__) || defined(__GNUC__)
#	define MND_UNSAFE_FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#	define MND_UNSAFE_FORCE_INLINE __forceinline
#else
#	define MND_UNSAFE_FORCE_INLINE inline
#endif

MND_UNSAFE_FORCE_INLINE constexpr f64 make_f64(u64 payload) noexcept {
	return __builtin_bit_cast(f64, payload);
}
MND_UNSAFE_FORCE_INLINE constexpr f32 make_f32(u32 payload) noexcept {
	return __builtin_bit_cast(f32, payload);
}

#undef MND_UNSAFE_FORCE_INLINE
} // namespace unsafe

#define MND_OPTION_ADD_CHECKED_ACCESS  \
	inline constexpr decltype(auto) value() & { \
		if(!has_value()) throw std::bad_optional_access(); \
		return (val_); \
	} \
	inline constexpr decltype(auto) value() const& { \
		if(!has_value()) throw std::bad_optional_access(); \
		return (val_); \
	} \
	inline constexpr decltype(auto) value() && { \
		if(!has_value()) throw std::bad_optional_access(); \
		return std::move(val_); \
	} \
	inline constexpr decltype(auto) value() const && { \
		if(!has_value()) throw std::bad_optional_access(); \
		return std::move(val_); \
	} \

#define MND_OPTION_ADD_UNCHECKED_ACCESS \
    inline constexpr decltype(auto) operator*() & \
		noexcept { return (val_); } \
    inline constexpr decltype(auto) operator*() const& \
		noexcept { return (val_); } \
    inline constexpr decltype(auto) operator*() && \
		noexcept { return std::move(val_); } \
    inline constexpr decltype(auto) operator*() const&& \
		noexcept { return std::move(val_); } \

inline constexpr u64 none_bits_f64 = 0x7ff80BADC0DEBEEFULL;
inline constexpr u32 none_bits_f32 = 0x7fcBADAFU;

inline constexpr f32 none_f32 = unsafe::make_f32(none_bits_f32);
inline constexpr f64 none_f64 = unsafe::make_f64(none_bits_f64);

/* We encode the enum discriminant in one of the quiet NANs.
 * Here we implicitly *trust* the users not to directly encode
 * std::nan(?) into the value, as this violates the discriminant.
 * For efficiency sake, this check is not enforced at runtime. */
#define IMPL_NICHE_FLOAT(N) \
struct pf##N { \
	using value_type = f##N; \
	\
	constexpr pf##N() : val_(none_f##N) {} \
	constexpr pf##N(None_t ) : val_(none_f##N) {} \
	constexpr pf##N(value_type val) : val_(val) {} \
	\
	inline constexpr bool has_value() const noexcept { \
		return __builtin_bit_cast(u##N, val_) != none_bits_f##N; \
	} \
	inline constexpr void reset() noexcept { val_ = none_f##N; } \
	\
	inline constexpr decltype(auto) value() & { \
		if(!has_value()) throw std::bad_optional_access(); \
		return (val_); \
	} \
	inline constexpr decltype(auto) value() const& { \
		if(!has_value()) throw std::bad_optional_access(); \
		return (val_); \
	} \
	inline constexpr decltype(auto) value() && { \
		if(!has_value()) throw std::bad_optional_access(); \
		return std::move(val_); \
	} \
	inline constexpr decltype(auto) value() const && { \
		if(!has_value()) throw std::bad_optional_access(); \
		return std::move(val_); \
	} \
	inline constexpr decltype(auto) operator*() &       noexcept { return (val_); } \
	inline constexpr decltype(auto) operator*() const&  noexcept { return (val_); } \
	inline constexpr decltype(auto) operator*() &&      noexcept { return std::move(val_); } \
	inline constexpr decltype(auto) operator*() const&& noexcept { return std::move(val_); } \
	\
private: \
	f##N val_; \
};

IMPL_NICHE_FLOAT(64)
IMPL_NICHE_FLOAT(32)
#undef IMPL_NICHE_FLOAT

/* Works for all signed/unsigned types, but also for non-owning pointers.
 * Now, shared pointer stuff I rather not touch, and honestly do not want to.
 * It is a can of worms, if put behind smart pointers, how to tell the smart pointer
 * which dtor it needs to invoke.. */
template<typename T>
struct NonZero {
	using value_type = T;

	constexpr NonZero() : val_(nval_) {};
	constexpr NonZero(None_t ) : val_(nval_) {};
	constexpr NonZero(T val) : val_(val) {};
	
	constexpr bool has_value() const noexcept { return val_ != nval_; }
	constexpr void reset() noexcept { val_ = nval_; }
	
	MND_OPTION_ADD_CHECKED_ACCESS
	MND_OPTION_ADD_UNCHECKED_ACCESS

private:
	T val_;
	constexpr static T nval_ = static_cast<T>(0);
};

} // namespace detail

using NonZeroU16 = detail::NonZero<u16>;
using NonZeroU32 = detail::NonZero<u32>;
using NonZeroU64 = detail::NonZero<u64>;
using NonZeroI16 = detail::NonZero<i16>;
using NonZeroI32 = detail::NonZero<i32>;
using NonZeroI64 = detail::NonZero<i64>;

template<typename T>
class Option : public detail::OptionBase<Option<T>, std::optional<T>> {
	using base = detail::OptionBase<Option<T>, std::optional<T>>;
public:
	using base::base;
	using base::operator=;
};

template<>
class Option<f64> : public detail::OptionBase<Option<f64>, detail::pf64> {
	using base = detail::OptionBase<Option<f64>, detail::pf64>;
public:
	using base::base;
	using base::operator=;
};

template<>
class Option<f32> : public detail::OptionBase<Option<f32>, detail::pf32> {
	using base = detail::OptionBase<Option<f32>, detail::pf32>;
public:
	using base::base;
	using base::operator=;
};

template<typename T>
class Option<T*> : public detail::OptionBase<Option<T*>, detail::NonZero<T*>> {
	using base = detail::OptionBase<Option<T*>, detail::NonZero<T*>>;
public:
	using base::base;
	using base::operator=;
};

template<typename T>
class Option<detail::NonZero<T>>
	: public detail::OptionBase<Option<detail::NonZero<T>>, detail::NonZero<T>> {
	using base = detail::OptionBase<Option<detail::NonZero<T>>, detail::NonZero<T>>;
public:
	using base::base;
	using base::operator=;
};

} // namespace mnd
