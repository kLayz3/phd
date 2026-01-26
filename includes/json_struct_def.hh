/* Few preproc black magic macros to help parse a JSON into a proper C++ structure.
 * It only works really for single-depth JSON where each entry's value is either
 * a string, integer or double, or array of these.
 * If the value is another JSON table, then can just encode it as a raw string.
 * We enhance the structure with some sprinkles to define structured binding.
 */

#pragma once

#include <iostream>

#define EMPTY_MACRO__(...) 

#define GET_HELP_AUX_IMPL  \
	template<std::size_t I> \
    decltype(auto) get() &       noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() const & noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() &&      noexcept { return get_helper<I>(std::move(*this)); } \
    template<std::size_t I> \
    decltype(auto) get() const&& noexcept { return get_helper<I>(std::move(*this)); } \

#define ADD_SERIALIZABLE_FIELD(TYPE, NAME, DEFAULT, INDEX) \
	TYPE NAME = DEFAULT; \
	using type_##INDEX = TYPE; \
	static constexpr const char* name_##INDEX = #NAME; \
	\
	template<std::size_t I, typename Self, typename std::enable_if<(I == INDEX)>::type* = nullptr> \
	static decltype(auto) get_helper(Self&& self) noexcept { return ( std::forward<Self>(self).NAME ); } \
	\
	template<std::size_t I, typename Self, typename std::enable_if<(I == INDEX)>::type* = nullptr> \
	static void print_field(std::ostream& os, Self&& self) noexcept { os << ".\e[1;94m" << #NAME << "\e[0m = " << std::forward<Self>(self).NAME << ", "; };

#define UNROLL_JSON_PARAM_SINGLE_(StructInstance, JSONInstance, INDEX) \
	try { \
		StructInstance.get<INDEX>() = JSONInstance[ decltype(StructInstance)::name_##INDEX ] \
			.get< \
				std::remove_reference_t< \
					decltype(StructInstance)::type_##INDEX \
				> \
			>(); \
	} catch(std::exception const& e) { \
		ERROR("Failed setup assignment \'%s\': index: %d, key:%s: reason: %s\n",  \
			#StructInstance, INDEX, decltype(StructInstance)::name_##INDEX, e.what()); \
	} \
	EMPTY_MACRO__(INDEX) \

#define UNROLL_JSON_PARAM(StructInstance, JSONInstance, N)  UNROLL_JSON_PARAM_N_(StructInstance, JSONInstance, N)
#define UNROLL_JSON_PARAM_N_(StructInstance, JSONInstance, N) UNROLL_JSON_PARAM_##N(StructInstance, JSONInstance)

#define UNROLL_JSON_PARAM_0(T1, T2)  EMPTY_MACRO__(T1, T2)     UNROLL_JSON_PARAM_SINGLE_(T1, T2, 0)
#define UNROLL_JSON_PARAM_1(T1, T2)  UNROLL_JSON_PARAM_0(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 1)
#define UNROLL_JSON_PARAM_2(T1, T2)  UNROLL_JSON_PARAM_1(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 2)
#define UNROLL_JSON_PARAM_3(T1, T2)  UNROLL_JSON_PARAM_2(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 3)
#define UNROLL_JSON_PARAM_4(T1, T2)  UNROLL_JSON_PARAM_3(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 4)
#define UNROLL_JSON_PARAM_5(T1, T2)  UNROLL_JSON_PARAM_4(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 5)
#define UNROLL_JSON_PARAM_6(T1, T2)  UNROLL_JSON_PARAM_5(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 6)
#define UNROLL_JSON_PARAM_7(T1, T2)  UNROLL_JSON_PARAM_6(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 7)
#define UNROLL_JSON_PARAM_8(T1, T2)  UNROLL_JSON_PARAM_7(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 8)
#define UNROLL_JSON_PARAM_9(T1, T2)  UNROLL_JSON_PARAM_8(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 9)
#define UNROLL_JSON_PARAM_10(T1, T2) UNROLL_JSON_PARAM_9(T1, T2) UNROLL_JSON_PARAM_SINGLE_(T1, T2, 10)

#define ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, INDEX) \
	namespace std { \
		template<> struct tuple_element<INDEX, TYPE> { using type = typename TYPE::type_##INDEX; }; \
	} \

#define ADD_STD_TYPE_RESOLUTION_(TYPE, N) \
	/* First define printer. */ \
	template<std::size_t... Is> \
	void print_me_json_helper_aux_##TYPE(std::ostream& os, const TYPE& obj, std::index_sequence<Is...>) { \
		(..., obj.print_field<Is>(os, obj) ); \
	} \
	inline std::ostream& operator<<(std::ostream& os, const TYPE& obj) { \
		os << "{"; \
		print_me_json_helper_aux_##TYPE(os, obj, std::make_index_sequence<N+1>{}); \
		os << "\e[1mEND\e[0m }"; \
		return os; \
	} \
	\
	namespace std { \
		template<> struct tuple_size<TYPE> : integral_constant<size_t, N+1> {};	\
	} \
	ADD_STD_TYPE_RESOLUTION_N_(TYPE, N) \

#define ADD_STD_TYPE_RESOLUTION_N_(TYPE, N) ADD_STD_TYPE_RESOLUTION_##N(TYPE)
#define ADD_STD_TYPE_RESOLUTION_0(TYPE)  EMPTY_MACRO__(TYPE)              ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 0)
#define ADD_STD_TYPE_RESOLUTION_1(TYPE)  ADD_STD_TYPE_RESOLUTION_0(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 1)
#define ADD_STD_TYPE_RESOLUTION_2(TYPE)  ADD_STD_TYPE_RESOLUTION_1(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 2)
#define ADD_STD_TYPE_RESOLUTION_3(TYPE)  ADD_STD_TYPE_RESOLUTION_2(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 3)
#define ADD_STD_TYPE_RESOLUTION_4(TYPE)  ADD_STD_TYPE_RESOLUTION_3(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 4)
#define ADD_STD_TYPE_RESOLUTION_5(TYPE)  ADD_STD_TYPE_RESOLUTION_4(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 5)
#define ADD_STD_TYPE_RESOLUTION_6(TYPE)  ADD_STD_TYPE_RESOLUTION_5(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 6)
#define ADD_STD_TYPE_RESOLUTION_7(TYPE)  ADD_STD_TYPE_RESOLUTION_6(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 7)
#define ADD_STD_TYPE_RESOLUTION_8(TYPE)  ADD_STD_TYPE_RESOLUTION_7(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 8)
#define ADD_STD_TYPE_RESOLUTION_9(TYPE)  ADD_STD_TYPE_RESOLUTION_8(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 9)
#define ADD_STD_TYPE_RESOLUTION_10(TYPE) ADD_STD_TYPE_RESOLUTION_9(TYPE)  ADD_STD_TYPE_RESOLUTION_SINGLE_(TYPE, 10)

