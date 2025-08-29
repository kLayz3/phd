#pragma once

#include "TOnceBase.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "TFile.h"
#include "libs.hh"
#include <cstring>
#include <type_traits>

class TH1;

namespace util {
	template<typename, typename = std::void_t<>>
	struct has_value_type : std::false_type {};

	template<typename U>
	struct has_value_type<U, std::void_t<typename U::value_type>> : std::true_type {};

	template<typename T>
	struct is_std_array : std::false_type{};

	template<typename T, std::size_t N>
	struct is_std_array<std::array<T, N>> : std::true_type {};

	template <typename T, typename = void>
	struct has_clone : std::false_type {};

	template <typename T>
	struct has_clone<T, std::void_t<decltype(std::declval<const T&>().Clone())>> : std::true_type {};

	template <typename T, typename = void>
	struct has_setname : std::false_type {};
 
	template <typename T>
	struct has_setname<T, std::void_t<
		decltype(std::declval<const T&>().SetName(std::declval<const char*>()))
		>> : std::true_type {};
}

/* A wrapper type for objects (such as TVectorD, THXX, TArray, TCutG, etc)
 * but also stl-containers (std::vector<int>, std::array<double, T>, etc) that might carry a name or maybe not.
 * This class expands them so that all of these classes' instances carry a name,
 * which shall be serialized (once) into a `TFile`,
 * unlike `TContainer`'s other parts which will be serialized row-wise into a `TTree` inside `TFile`. */

template<typename T>
class TOnce : public TOnceBase {
    static_assert(! std::is_pointer_v<T>, "Must not pass pointer type (T*) to TOnce<T>");
    static_assert(! std::is_fundamental_v<T>, "Must not pass trivial type. Wrap it in e.g. TParameter<T> first.");
    static_assert(! std::is_void_v<T>, "Hello?");
    static_assert(! std::is_array_v<T>, "Must not pass raw C-style arrays. Pass an `std::array<T,N>` instead.");
	static_assert(  std::is_copy_assignable_v<T>, "Type T must be copy-assignable");

    T _internal;

public:
    using type = T;

    /* Case 1: T constructible with (const char*, Args...) */
    template<typename... Ts,
        typename std::enable_if<std::is_constructible_v<T, const char*, Ts...>>::type* = nullptr
            > TOnce(const char* name, Ts&&... args) : TOnceBase(name), 
            _internal(name, std::forward<Ts>(args)...) {}

    /* Case 2: T constructible with (Args...) but not (const char*, Args...) */
    template<typename... Ts,
        typename std::enable_if<
            !std::is_constructible_v<T, const char*, Ts...> &&
            std::is_constructible_v<T, Ts...>
            >::type* = nullptr
        > TOnce(const char* name, Ts&&... args) : TOnceBase(name),
            _internal(std::forward<Ts>(args)...) {}

    /* Weird case: T constructible via init list; C++ painpoint. */
    template<typename U = T, typename std::enable_if<util::has_value_type<U>::value>::type* = nullptr>
        TOnce(const char* name, std::initializer_list<typename U::value_type> il) :
            TOnceBase(name) {
				if constexpr(std::is_constructible_v<U, std::initializer_list<typename U::value_type>>)
					_internal = U(il); // vector, list, etc
				else if constexpr(util::is_std_array<U>::value) {
					if(il.size() == 0)
						_internal.fill(typename U::value_type() );
					else if(il.size() == 1)
						_internal.fill(*il.begin());
					else if(il.size() == _internal.size())
						std::copy(il.begin(), il.end(), _internal.begin());
					else
						assert(il.size() == _internal.size() && "Initializer list for std::array<T,N> must have either 0, 1 or N members"); 
				}
				else static_assert(std::is_constructible_v<U, std::initializer_list<typename U::value_type>>,
						"Type doesn't support initializer lists ctor.");
			}

    template<typename std::is_default_constructible<T>::type* = nullptr>
        TOnce(const char* name = "") : TOnceBase(name), _internal() {}
	
	TOnce(const TOnce& ) = default;
	TOnce& operator=(const TOnce& ) = default;
	TOnce(TOnce&& ) noexcept = default;
	TOnce& operator=(TOnce&& ) noexcept = default;
	~TOnce() = default;

	void SetName(std::string name, const char* title = "") {
        TOnceBase::SetName(name);
		if constexpr(std::is_base_of_v<TNamed, T> || util::has_setname<T>::value) {
			_internal.SetName(name.c_str());
			if constexpr(std::is_base_of_v<TNamed, T>)
				if(title && *title) _internal.SetTitle(title);
		}
    }

    Int_t Write(TFile* f = nullptr, const char* name = "") override {
		if(!name) ERROR("Don't pass nullptr for `name` here.");
		if(!f || f->IsZombie() || !f->IsOpen()) {
			f = gDirectory->GetFile();
			if(!f || f->IsZombie() || !f->IsOpen())
				ERROR("%s (label: \'%s\') : output ROOT file not supplied or invalid and gDirectory holds no open valid file.", _name.c_str(), name);
		}
        if constexpr(std::is_base_of_v<TObject, T>)
            return f->WriteTObject(&_internal, *name ? name : _name.c_str());
        else 
            return f->WriteObject (&_internal, *name ? name : _name.c_str());
    }
    
	void* Load(TFile* f = nullptr, const char* target = "") override {
		const char* name = *target ? target : _name.c_str();
		if(!name || ! *name) ERROR("Unnamed TOnce<T> object while trying to load from a file.");
		if(!f || f->IsZombie() || !f->IsOpen()) {
			f = gDirectory->GetFile();
			if(!f || f->IsZombie() || !f->IsOpen()) 
				ERROR("(%s) (target: \'%s\') tried to open the file, but didn't receive a handle and also gDirectory holds no open valid file.", _name.c_str(), target);
		}
		T* tmp = f->Get<T>(name);
		if(!tmp) ERROR("(%s) - asked for " EMPH(%s) " name as key to a static object, got back nullptr.", _name.c_str(), name);
		try {
			_internal = *tmp; // deep-copy ; asserted via type traits on the top.
			if constexpr(std::is_base_of_v<TH1, T>) { 
				tmp->SetDirectory(nullptr); // remove from ROOT’s ownership
				delete tmp;
			}
		} catch(const std::exception& e) {
			ERROR("Caught: " EMPH(%s) " while trying to copy during loading of a static object.", e.what());	
		} catch(...) {
			ERROR("Unknown exception caught.");
		}
		return (void*)&_internal;
	}

	inline T& operator()() { return _internal; }
	inline const T& operator()() const { return _internal; }

	inline T* operator->() { return &_internal; }
	inline const T* operator->() const { return &_internal; }
};

