#pragma once

#include <unordered_map>
#include "AuxFunctions.hh"

using TDictInfo = std::unordered_map<std::string, std::string>;

struct TContainerBase {
	std::string _name;
	
	inline const char* GetName() const noexcept { return _name.c_str(); }
	inline void SetName(std::string(name)) noexcept { _name = std::move(name); }
	
	TContainerBase(std::string name) : _name(std::move(name)) {}

	/** 
	 * Should be initial (optional) call. Setting up metadata, names and delegating ctors to objects.
	 * Override it in the derived classes.
	 */
	inline void Init(TDictInfo info) { (void)info; }
};

/**
 * Encapsulates a type needed to be used as a (de)serialization target
 * from/to RNTuple column. 
 */
template<typename T>
struct TContainer : TContainerBase {
	static_assert( std::is_copy_constructible_v<T>, "Type must be copy constructible.");
	static_assert( std::is_copy_assignable_v<T>, "Type must be copy assignable.");
	static_assert( std::is_move_constructible_v<T>, "Type must be move constructible.");
	static_assert( std::is_move_assignable_v<T>, "Type must be move assignable.");
	static_assert(!std::is_void_v<T>, "Hello?");
	static_assert(!std::is_array_v<T>, "Must not pass raw C-style arrays. Pass an `std::array<T,N>` instead.");
	static_assert(!std::is_pointer_v<T>, "Must not pass pointer type (T*).");
	static_assert(!std::is_reference_v<T>, "Must not pass ref type (T&).");

	using inner_type = T;

private:
	std::shared_ptr<T> _inner;

public:
	TContainer() {};
	TContainer(std::string name) : TContainerBase(name) {}
	
	T& operator*() noexcept { return *(_inner.get()); }
	T& inner() noexcept	{ return this->operator*(); }

	void Clean() noexcept {
		if constexpr(util::has_clean_noexcept<T>::value)
			_inner->Clean();
		else
			static_assert(util::has_clean_noexcept<T>::value, "Type has no `void Clean() noexcept` method.\n");
	}
};

/* Type `T` is either something like TXXXYYYEvent (Go4) or a custom structure
 * for UCESB. 
 * Namely, UCESB splits all the trees directly into leaves, so address has to be mapped
 * sequentally.
 * Encapsulates a type needed to be used as a deserialization target from a TTree branch. 
 */
template<typename T>
struct TRawContainer : TContainerBase {
	static_assert(!std::is_void_v<T>, "Hello?");
	static_assert(!util::is_an_array_v<T>, "Must not pass raw C-style array or `std::array<T,N>`.");
	static_assert(!std::is_pointer_v<T>, "Must not pass pointer type (T*).");
	static_assert(!std::is_reference_v<T>, "Must not pass ref type (T&).");
	
	using inner_type = T;

private:
	T* _inner;

public:
	TRawContainer() = default;
	TRawContainer(std::string name) : TContainerBase(name) {}

	T& operator*() noexcept { return *_inner; }
	T& inner() noexcept	{ return this->operator*(); }

};
