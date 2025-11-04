#pragma once

#include <chrono>
#include <cstdio>
#include <sstream>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <variant>
#include <vector>
#include "libs.hh"
#include <fstream>
#include <regex>
#include <typeinfo>
#include <optional>

#include "nlohmann/json.hpp"

#ifndef _MSC_VER
#	include <cxxabi.h>
#endif

class TObject;

#if __has_include("indicators/indicators.hh")
#	define __HAS_INDICATORS
#	include "indicators/indicators.hh"
#endif

namespace util {

#ifdef __HAS_INDICATORS
inline void 
PrintProgress(indicators::ProgressBar& bar,	const u64 n_entry, const u64 max_entries, const u64 step = 250) noexcept {
	static u64 n_entry_called = 0;
	if(n_entry - n_entry_called < step) return;
	bar.set_progress( (n_entry*100) / max_entries );

	n_entry_called = n_entry;
}
#endif

template<typename T>  
void QuickSwap(std::vector<T>& v, int i, int j) noexcept {
	if(!v.size()) return;
    std::swap(v[i], v[j]);
}

template<typename T>  
void QuickErase(std::vector<T> &v, int i) noexcept {
    std::swap(v[i], v.back());
    v.pop_back();
}

template<typename T>
void Append(std::vector<T>& dst, const std::vector<T>& src) noexcept {
	dst.reserve(dst.size() + src.size());
	dst.insert(dst.end(), src.begin(), src.end());
}

template<typename T>
void Append(std::vector<T>& dst, std::vector<T>&& src) noexcept {
	if(&dst == &src) return;
	dst.reserve(dst.size() + src.size());
	dst.insert(dst.end(),
			std::make_move_iterator(src.begin()),
			std::make_move_iterator(src.end()));
	/* src vector left in undefined state. */
}

/* ------------------------- */
template<typename T>
T rround(double x) noexcept { return static_cast<T>(x + 0.5); }

using std::chrono::duration_cast;
using std::chrono::seconds;
using std::chrono::milliseconds;

#define timeNow() std::chrono::high_resolution_clock::now()
struct TimePoint {
	std::chrono::high_resolution_clock::time_point t;
	std::string tag;
	TimePoint() : t(timeNow()), tag("?") {};
	TimePoint(std::string s) : t(timeNow()), tag(s) {};
};
#undef timeNow

enum TimingVariant : u64 {
	kMINUTE = 1,
	kSECOND = 60,
	kMILLISECOND = 60'000,
	kMICROSECOND = 60'000'000,
};

template<TimingVariant E = kMILLISECOND>
void PrintElapsed(const TimePoint& end, const TimePoint& start) {
	using us = std::chrono::microseconds;

	static const char* mode[] = {"us", "ms", "sec", "min"};
	int mode_i = 0;
	double elapsed;
	if constexpr(E == kMINUTE) {
		elapsed = (std::chrono::duration_cast<us>(end.t - start.t).count()) * static_cast<double>(kMINUTE) / kMICROSECOND;
		mode_i = 3;
	}
	else if constexpr(E == kSECOND) {
		elapsed = (std::chrono::duration_cast<us>(end.t - start.t).count()) * static_cast<double>(kSECOND) / kMICROSECOND;
		if(elapsed > 60) {
			elapsed *= static_cast<double>(kMINUTE) / kSECOND;
			mode_i = 3;
		} else 
			mode_i = 2;
	}
	else if constexpr(E == kMILLISECOND) {
		elapsed = (std::chrono::duration_cast<us>(end.t - start.t).count()) * static_cast<double>(kMILLISECOND) / kMICROSECOND;
		mode_i = 1;
	}
	else if constexpr(E == kMICROSECOND) {
		elapsed = (std::chrono::duration_cast<us>(end.t - start.t).count());
		mode_i = 0;
	}
    printf("Elapsed time from " EMPH1(%s) " to " EMPH1(%s) ": " EMPH(%.1f) " %s\n", start.tag.c_str(), end.tag.c_str(), elapsed, mode[mode_i]);
}

template<TimingVariant E = kMILLISECOND>
inline void PrintElapsed(const std::vector<TimePoint>& v) {
	if(v.size() < 2) return;
	PrintElapsed<E>(v[v.size()-1], v[v.size() - 2]);
}

template<TimingVariant E = kMILLISECOND>
inline void PrintElapsed(std::vector<TimePoint>&& v) {
	if(v.size() < 2) return;
	printf("Total execution time: ");
	PrintElapsed<E>(v.back(), v.front());
}


inline void __concat_impl__(std::ostringstream& ) {}

template<typename T, typename... Args>
void __concat_impl__(std::ostringstream& oss, T&& first, Args&&... args) {
	oss << std::forward<T>(first);
	__concat_impl__(oss, std::forward<Args>(args)... );
}

/**
 * Concatenates bunch of arguments which can be lvalues, statics, etc. and returns an owned std::string. 
 */
template<typename... Args>
std::string sstrcat(Args&&... args) {
	std::ostringstream oss;
	__concat_impl__(oss, std::forward<Args>(args)... );
	return oss.str();
}

/* ----- std::variant trickery. ----- */
using Empty = std::monostate;

template<typename T>
using Maybe = std::optional<T>;

template<typename... Ts>
using Variant = std::variant<Empty, Ts...>;

template<typename... Ts>
constexpr bool IsEmpty(const Variant<Ts...>& v) {
	return v.valueless_by_exception() or std::holds_alternative<Empty>(v);
}

/**
 * Apply 'std::visit' https://en.cppreference.com/w/cpp/utility/variant/visit2.html
 * to the util::Variant type. Throws if the 'util::Empty' type is encountered. 
 */
template<typename Visitor, typename Var>
decltype(auto) visit_non_empty(Visitor&& vs, Var&& /* Variant&& */ var) {
	return std::visit(
		[&](auto&& value) -> decltype(auto) {
			using T = std::decay_t<decltype(value)>;
			if constexpr(std::is_same_v<T, Empty>)
				throw std::runtime_error("util::visit_non_empty except: `Empty` type encountered.");
			else 
				return std::invoke(std::forward<Visitor>(vs),
					std::forward<decltype(value)>(value));
		},
		std::forward<Var>(var)
	);
}

template<typename Tuple, typename Callable, std::size_t... Is>
void _for_each_in_tuple_impl(Tuple&& t, Callable&& f, std::index_sequence<Is...>) {
	(..., f (std::get<Is>(std::forward<Tuple>(t))));
}

template<typename Tuple, typename Callable>
void for_each_in_tuple(Tuple&& t, Callable&& f) {
	constexpr std::size_t N = std::tuple_size_v<std::decay_t<Tuple>>;
	_for_each_in_tuple_impl(std::forward<Tuple>(t), std::forward<Callable>(f), 
		std::make_index_sequence<N>{});
}

template<typename T>
struct is_an_array : std::false_type {};
template<typename T, size_t N>
struct is_an_array<T[N]> : std::true_type { using value_type = T; static constexpr size_t size = N; };
template<typename T, size_t N>
struct is_an_array<std::array<T,N>> : std::true_type { using value_type = T; static constexpr size_t size = N; };
template<typename T>
inline constexpr bool is_an_array_v = is_an_array<T>::value;

template<typename T, typename = void>
struct has_clean_noexcept : std::false_type {};
template<typename T>
struct has_clean_noexcept<T, 
	std::void_t<decltype(std::declval<T&>().Clean())>>
	: std::bool_constant<
		std::is_same_v<void, decltype(std::declval<T&>().Clean())> &&
		noexcept(std::declval<T&>().Clean())
	> {};

template<typename, typename = std::void_t<>>
struct has_value_type : std::false_type {};
template<typename U>
struct has_value_type<U, std::void_t<typename U::value_type>> : std::true_type {};

template<typename T>
struct is_std_array : std::false_type{};
template<typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template<typename T, typename = void>
struct has_clone : std::false_type {};
template<typename T>
struct has_clone<T, std::void_t<decltype(std::declval<const T&>().Clone())>> : std::true_type {};

template<typename T, typename = void>
struct has_setname : std::false_type {};
template<typename T>
struct has_setname<T, std::void_t<
	decltype( std::declval<T&>().SetName(std::declval<const char*>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_copy : std::false_type {};
template<typename T>
struct has_copy<T, std::void_t<
	decltype( std::declval<const T&>().Copy( std::declval<TObject&>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_set_directory : std::false_type {};
template<typename T>
struct has_set_directory<T, std::void_t<
	decltype( std::declval<const T&>().SetDirectory( std::declval<const char*>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_dyadic_add_ref : std::false_type {};
template<typename T>
struct has_dyadic_add_ref<T, std::void_t<
	decltype( std::declval<T&>().Add( std::declval<const T&>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_dyadic_add_ptr : std::false_type {};
template<typename T>
struct has_dyadic_add_ptr<T, std::void_t<
	decltype( std::declval<T&>().Add( std::declval<const T*>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_free_add_fn : std::false_type {};
template<typename T>
struct has_free_add_fn<T, std::void_t<
	decltype( Add( std::declval<T&>(), std::declval<const T&>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_free_mean_fn : std::false_type {};
template<typename T>
struct has_free_mean_fn<T, std::void_t<
	decltype( Mean( std::declval<T&>(), std::declval<const T&>() ) )
>> : std::true_type {};

template<typename T, typename = void>
struct has_add_div_unary_op : std::false_type {};
template<typename T>
struct has_add_div_unary_op<T, std::void_t <
	std::conjunction <
		decltype( std::declval<T&>().operator+=(std::declval<const T&>() ) ) ,
		decltype( std::declval<T&>().operator/=(std::declval<const int>() ) )
	>
>> : std::true_type {};

template< template<typename...> typename Base, typename Derived>
struct _is_base_of_template_impl {
	template<typename... Ts>
	static constexpr std::true_type test(Base<Ts...>& );
	static constexpr std::false_type test(...);
	using type = decltype(test(std::declval<Derived&>()));
};

template< template<typename...> typename Base, typename Derived>
using is_base_of_template = typename _is_base_of_template_impl<Base, Derived>::type;

template<typename Tuple, std::size_t... Is>
auto zip_refs_impl(Tuple& t1, const Tuple t2, std::index_sequence<Is...>) {
	using R = std::tuple <
		std::pair <
			std::add_lvalue_reference_t<std::tuple_element_t<Is, Tuple>>,
			std::add_lvalue_reference_t<std::add_const_t<std::tuple_element_t<Is, Tuple>>>
		>...
	>;
	
	return R {
		std::pair <
			std::add_lvalue_reference_t<std::tuple_element_t<Is, Tuple>>,
			std::add_lvalue_reference_t<std::add_const_t<std::tuple_element_t<Is, Tuple>>>
		> { std::get<Is>(t1), std::get<Is>(t2) }...
	};
}

template<typename Tuple>
auto zip_refs(Tuple& t1, const Tuple& t2) {
	constexpr std::size_t N = std::tuple_size_v<Tuple>;
	return zip_refs_impl( t1, t2, std::make_index_sequence<N>{} );
}

template<typename T>
constexpr bool is_pathlike_arg_v =
	std::is_base_of_v<std::ifstream, std::decay_t<T>> ||
	std::is_constructible_v<std::ifstream, T&&> ||
	std::is_constructible_v<std::string, T&&> ||
	std::is_constructible_v<std::filesystem::path, T&&>;

template<typename T,
	typename U = std::decay_t<T>,
	typename std::enable_if<is_pathlike_arg_v<T>>::type* = nullptr
> std::optional<std::ifstream> get_maybe_ifstream(T&& arg) {
	if constexpr(std::is_base_of_v<std::istream, U>) {
		if(arg.is_open()) return arg;
		else return {};
	}
	else if constexpr(
		std::is_constructible_v<std::ifstream, T&&> ||
		std::is_constructible_v<std::filesystem::path, T&&>
	) {
		auto f = std::ifstream(std::forward<T>(arg));
		if(f.is_open()) return f;
		else return {};
	}
	else if constexpr(std::is_constructible_v<std::string, T&&>) {
		auto f = std::ifstream(std::string(std::forward<T>(arg)));
		if(f.is_open()) return f;
		else return {};
	}
	else 
		static_assert(std::is_constructible_v<std::string, T&&>,
			"Type not constructible to ifstream or path-like type.");
}

template<typename T>	
bool is_file_readable(T&& arg) {
	return get_maybe_ifstream(std::forward<T>(arg)).has_value();
}

template<typename T,
	typename U = std::decay_t<T>
> std::optional<std::string> get_file_path(T&& arg) {
	static_assert(is_pathlike_arg_v<U>, "Type <T> must be path-like");
	if constexpr(std::is_base_of_v<std::istream, U>) 
		return {};
	else if constexpr(std::is_constructible_v<std::string, T&&>)
		return std::string(std::forward<T>(arg));
	else if constexpr(std::is_constructible_v<std::filesystem::path, T&&>)
		return std::filesystem::path(std::forward<T>(arg)).string();
}

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>,
	std::size_t N = is_an_array<U>::size,
	char(*)[N % 2] = nullptr // Odd
> auto median(const T& sorted_arr) { return sorted_arr[(N-1)/2]; }

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>,
	std::size_t N = is_an_array<U>::size,
	char(*)[!(N % 2)] = nullptr // Even
> auto median(const T& sorted_arr) { return ( sorted_arr[(N-1)/2] + sorted_arr[N/2] ) / 2; }

template<typename T, typename... Ts>
constexpr auto min(T t, Ts... ts) noexcept {
	static_assert(std::is_arithmetic_v<
		std::remove_reference_t<T>
	> &&
	(std::is_arithmetic_v<
		std::remove_reference_t<Ts>
	> && ...), "Types passed must be arithmetic, basic types.");

	using R = std::common_type_t<T, Ts...>;
	R r = static_cast<R>(t);
	((r = std::min<R>(r, static_cast<R>(ts))), ...);
	return r;
}
template<typename T, typename... Ts>
constexpr auto max(T t, Ts... ts) noexcept {
	static_assert(std::is_arithmetic_v<
		std::remove_reference_t<T>
	> &&
	(std::is_arithmetic_v<
		std::remove_reference_t<Ts>
	> && ...), "Types passed must be arithmetic, basic types.");

	using R = std::common_type_t<T, Ts...>;
	R r = static_cast<R>(t);
	((r = std::max<R>(r, static_cast<R>(ts))), ...);
	return r;
}

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>
> constexpr int FindIndex(const T& arr, const typename is_an_array<U>::value_type& val) {
	static_assert(is_an_array_v<U>, "Type T must be a C-style array or \'std::array<T,N>\'");
	constexpr std::size_t N = is_an_array<U>::size;
	for(int i=0; i < (int)N; ++i)
		if(arr[i] == val) return i;
	return -1;
}

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>
> constexpr int len(const T& arr) {
	static_assert(is_an_array_v<U>, "Passed type must either be an array reference &(T[N]) or &std::array<T,N>.");
	return static_cast<int>(is_an_array<U>::size);	
}

static inline void Trim(std::string& s) {
	s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c){return std::isspace(c);}), s.end());
}

using nlohmann::json;

static inline void parse_json_string(std::vector<int>& out, std::string s) {
	static const std::regex re_num(
		R"(^(0x|0b)?(\d+)$)"
	);
	static const std::regex re_range(
		R"(^(\d+)\.\.(=)?(\d+)$)"
	);
	static const std::regex re_seq(
		R"(^(\d+)n(\+\d+)?$)"
	);

	Trim(s);

	/* Strings can be passed either as:
	 * 1) raw numbers
	 * 2) range x1..x2 (x2 excluded) or x1..=x2 (x2 included)
	 * 3) sequence: an+b ('a', 'b' are the parameters, '+b' optional; 'n' fixed token). 
	 *    Meaning: strips: b, a+b, 2*a+b, etc. */
	std::smatch m;
	if(std::regex_match(s, m, re_num)) {
		if(m[1].matched) {
			if(m[1].str() == "0x") out.push_back(std::stoi(m[2].str(), nullptr, 16));
			else out.push_back(std::stoi(m[2].str(), nullptr, 2));
		}
		else out.push_back(std::stoi(m[2].str()));
	} else if(std::regex_match(s, m, re_range)) {
		int a = std::stoi(m[1].str());
		int b = std::stoi(m[3].str());
		if(a > b) ERROR("%d < %d found while parsing json: \'%s\'\n", a,b,s.c_str());
		for(int i=a; i<b; ++i) out.push_back(i);
		if(m[2].matched) out.push_back(b);
	} else if(std::regex_match(s, m, re_seq)) {
		int a = std::stoi(m[1].str());
		int b = m[2].matched ? std::stoi(m[2].str()) : 0;
		for(int i=b; i <	 640; i+=a)
			out.push_back(i);
	} else {
		WARN("String \'%s\' doesn't match a: number, range or sequence regular expression.", s.c_str());
	}
}

/**
 * Json Custom range parser, accepting raw numbers or strings (or arrays of the same)
 * parsable as int or as a range:
 * 'a1..a2' left-inclusive, or 'a1..=a2' right inclusive.
 * */
static inline void parse_json_as_int_vec(std::vector<int>& out, const json& j) {
	if(j.is_string()) {
		parse_json_string(out, j.get<std::string>());
	}
	else if(j.is_number()) {
		out.push_back(j.get<int>());	
	}
	else if(j.is_array()) {
		out.reserve(j.size());
		for(const auto& jsub : j) 
			parse_json_as_int_vec(out, jsub);
	}
	else ERROR("Passed a json object '%s' which isn't: string/array/number.\n", j.dump().c_str());
}

static inline void append_flat_json(json& dst, const json& src) {
	if(!dst.empty() && !dst.is_object())
		ERROR("Destination object \'%s\' cannot store sources appended to it. Non-empty and not-an-object!", dst.dump().c_str());
	if(!src.is_object())
		ERROR("Source json: \'%s\' must be an object.", src.dump().c_str());

	for(const auto& [k, v] : src.items()) {
		auto it = dst.find(k);
		if(it == dst.end() || it->is_null()) {
			dst[k] = v;
			continue;
		}

		if(!it->is_array()) {
			*it = json::array({ *it }); // [old]
		}

		if(v.is_array()) {
			it->insert(it->end(), v.begin(), v.end()); // "k": [..., v[0], v[1], v[n-1] ]
		} else {
			it->push_back(v); // "k": [..., v]
		}
	}
}

#if __has_include(<new>)
#	include <new>
#endif

/* Returns the name of the type passed, also adding 
 * ref-cv qualifiers. Demangles templated types, too :-) */
template<typename T,
	typename U = std::decay_t<T>	
> std::string type_name() {
	std::unique_ptr<char, void(*)(void*)> own(
#ifndef _MSC_VER
		abi::__cxa_demangle(typeid(U).name(), nullptr,
			nullptr, nullptr),
#else
		nullptr,
#endif
		std::free
	);
	std::string r = (own) ? own.get() : typeid(U).name();
	if constexpr(std::is_const<T>::value)
        r += " const";
    if constexpr(std::is_volatile<T>::value)
        r += " volatile";
    if constexpr(std::is_lvalue_reference<T>::value)
        r += "&";
    else if constexpr(std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}

#define _SELF_TYPE_CSTR \
	util::type_name<typename std::remove_reference<decltype(*this)>::type>().c_str()

/* Sometimes cursor can be hidden mid execution,
 * if the program dies due to a system signal,
 * execute this to bring it back. Only POSIX async-safe 
 * calls are allowed. NOTE: this *might* disable standard ROOT
 * stack trace dump in case of SIGSEGV catch. Unmap it in this case.
 * and in Linux execute 'tput cnorm' to get the cursor back. */
inline void sig_callback_handler(int signum) {
	const char show[] = "\x1b[?25h";
    const char nl   = '\n';
	write(STDERR_FILENO, &nl, 1);
	write(STDERR_FILENO, show, sizeof show - 1);
	WARN_ASYNC("Caught abort/seg signal [%d].\n", signum);
	_exit(128 + signum);
}

constexpr std::size_t CL = 
#if defined(__cpp_lib_hardware_interference_size) && defined(__clang__)
	std::hardware_destructive_interference_size
#elif defined(__x86_64__)
	64 // fallback, for sure.
#else
	64 // Suppose.
#endif
	;

} // namespace util
