#pragma once

#include "Rtypes.h"
#include "libs.hh"
#include "AuxFunctions.hh"
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

template<typename T> class TContainer;

/**
 * Base class API for all types that will implement `void ProcessEntry()` and `Int_t Write()` methods.
 * It holds refs/ptrs to associated containers, or the containers should be movable.
 */

class TProcessorBase {
public:
	struct IOInfo { std::string file_name; std::string rntuple_name; };
	
	inline static IOInfo info_in {}, info_out {};

	TProcessorBase() = default;
	TProcessorBase(const TProcessorBase& ) = default;
	TProcessorBase& operator=(const TProcessorBase& other) = default;
	TProcessorBase(TProcessorBase&& other) = default;
	TProcessorBase& operator=(TProcessorBase&& other) = default;

	virtual ~TProcessorBase() = default;

public:
	inline virtual Int_t Write() {
		WARN("Write called but from the base class (no-op: Nothing is written). \n\tForgot to override it in derived class?\n");
		return 0; 
	}
};

struct TContainerBaseW;
struct TContainerBaseR;
template<typename T> struct TRawContainer;

template<typename>
struct TProcessor; /* Undefined. */

/**
 * Represents a single analysis subprocess, mapping the data
 * from combination of Input structures to one unique Output structure.
 */
template<typename Out, typename... Ins>
struct TProcessor<Out(Ins...)> : TProcessorBase {
	static_assert((std::disjunction_v<
			util::is_base_of_template<TContainer, Ins>,
			util::is_base_of_template<TRawContainer, Ins>> && ...), 
		"Input type(s) must inherit from (or be) TContainer<T> / TRawContainer<T>.");
	static_assert(util::is_base_of_template<TContainer, Out>::value, "Output type must inherit from (or be) TContainer<T>.");

	std::tuple<Ins...> in;
	Out out;
	
	/**
	 * Returns a reference to the input container.
	 */
	template<u32 N = 0> 
	decltype(auto) GetInput() const& {
		static_assert(N < std::tuple_size_v<decltype(in)>, "Accessing N >= tuple size.");
		return ( std::get<N>(in).inner() );
	}

	template<u32 N = 0> 
	decltype(auto) GetInput() & {
		static_assert(N < std::tuple_size_v<decltype(in)>, "Accessing N >= tuple size.");
		return ( std::get<N>(in).inner() );
	}
	
	TProcessor() = default;
	
	/** 
	 * Only other acceptable ctor is to take ownership from an existing output object.
	 * The output object need to be specifically `std::move`'ed to be perfectly clear in user code.
	 * Input objects, might be shared between different subprocesses - here we make our own copy.
	 */
	TProcessor(Out&& _out, const Ins&... ins) {
		this->out = std::move(_out);
		in = std::make_tuple(ins...);
	}

	static_assert(std::is_copy_constructible_v<TProcessor>,
		"TProcessor<Out(Ins...)> must be copy-constructible. are all the Container objects so?");
	static_assert(std::is_copy_assignable_v<TProcessor>,
		"TProcessor<Out(Ins...)> must be copy-assignable. are all the Container objects so?");
	static_assert(std::is_move_constructible_v<TProcessor>,
		"TProcessor<Out(Ins...)> must be move-constructible. are all the Container objects so?");
	static_assert(std::is_move_assignable_v<TProcessor>,
		"TProcessor<Out(Ins...)> must be move-assignable. are all the Container objects so?");
};
