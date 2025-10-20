#pragma once
#include "AuxFunctions.hh"

template<typename T>
struct TRawContainer {
	static_assert(!std::is_void_v<T>, "Hello?");
	static_assert(!std::is_array_v<T>, "Must not pass raw C-style arrays. Pass an `std::array<T,N>` instead.");
	static_assert(!std::is_pointer_v<T>, "Must not pass pointer type (T*).");
	static_assert(!std::is_reference_v<T>, "Must not pass ref type (T&).");

	std::string _name;
	T* _inner;

public:
	TRawContainer() = default;
	TRawContainer(std::string name) : TRawContainer() {
		_name = std::move(name);
	}

	T& operator*() noexcept { return *_inner; }
	T& inner() noexcept	{ return this->operator*(); }

	const char* GetName() const { return _name.c_str(); }
	void SetName(std::string(name)) { _name = std::move(name); }
};
