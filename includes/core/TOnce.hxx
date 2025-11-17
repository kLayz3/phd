#pragma once

#include "RtypesCore.h"
#include "TOnceBase.hxx"
#include "TDirectory.h"
#include "TNamed.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "libs.hh"
#include <cstring>
#include "AuxFunctions.hh"
#include <type_traits>

class TH1;
class TTree;
template<typename T> struct TContainer;

/**
 * A wrapper type for objects (such as TVectorD, THXX, TArray, TCutG, etc)
 * but also stl-containers such as std::vector<int>, std::array<double, T>, etc that might carry 
 * a name or maybe not. This class expands them so that all of these classes' instances carry a name,
 * which shall be serialized (once) into a `TFile`,
 * unlike `TContainer`'s other parts which will be serialized row-wise into the RNTuple.
 * Instances are differentiated by their unique string key `_name`.
 *
 * Available types to be wrapped are limited to only the ones that can be either Summed or Averaged together.
 * Or you supply a custom collector to the constructor of this wrapper type.
 */

/* Different threads will have their own unique exclusive objects of this wrapper type (for writing),
 * e.g. TH1I's. Once their work is done, there has to be a way to collect and combine them into a single object,
 * that is to be then written into the output ROOT file.
 * Here we differentiate two such categories of types.
 * [1] Sum type - T must implement one of the following, checked in order:
 *         void T::Add(const T& rhs)   (such as 'TH1')
 *         void Add(T& lhs, const T& rhs) (free function; custom types; e.g. std::array<T>)
 *         T& operator+=(const T& rhs) (such as 'int')
 * [2] Mean types - T must implement one of the following, checked in order:
 *         void Mean(T& lhs, const T& rhs)
 *         T& operator+=(const T& rhs) and T& operator/=(const int) (such as 'int', 'double')
 *
 *  The collector will sum up all the instances of identical Sum type, 
 *  and make a mean value of all the instances of identical Mean type.
 *  Collecting is done via a dyadic fold, e.g. for 8 instances of TOnce<T> objects:
 *
 *  a0, a1, a2, a3, a4, a5, a6, a7 
 *    \/      \/      \/      \/    
 *   a01,    a23,    a45,    a67
 *      \   /           \   /
 *      a0123     ,     a4567
 *           \         /
 *            a01234567
 *
 * It is done in-place, a0 now carries the summed/mean'ed up value, while the other (N-1)
 * instances carry possibly intermediate calculations, and should not be used. 
 * In the master Pool instantization, it's asserted that N forms a perfect log(2) (N == 2^n, for some n).
 */

template<typename T>
class TOnce : public TOnceBase {
	static_assert(! std::is_pointer_v<T>, "Mst not pass pointer type (T*) to TOnce<T>");
	static_assert(! std::is_fundamental_v<T>, "Must not pass trivial type. Wrap it in e.g. TParameter<T> first.");
	static_assert(! std::is_void_v<T>, "Hello?");
	static_assert(! std::is_array_v<T>, "Must not pass raw C-style arrays. Pass an `std::array<T,N>` instead.");
	static_assert(! std::is_base_of_v<TTree, T>, "Cannot wrap the TTree type here.");
	static_assert(  std::is_copy_assignable_v<T>, "Type T must be copy-assignable");
	template<typename U> friend struct TContainer;

	void (*_collector)(T&, const T&) = nullptr;
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
			!std::is_constructible_v<T, const char*, Ts...> && std::is_constructible_v<T, Ts...>
			>::type* = nullptr
		> TOnce(const char* name, Ts&&... args) : TOnceBase(name),
		_internal(std::forward<Ts>(args)...) {}

	/* Case 3: T constructible via init list; C++ painpoint. */
	template<typename U = T, typename std::enable_if<util::has_value_type<U>::value>::type* = nullptr>
		TOnce(const char* name, std::initializer_list<typename U::value_type> il) : TOnceBase(name) {
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
				"Type doesn't correctly support initializer lists ctor.");
		}

	/* Case 4: default constructor. */
	template<typename std::is_default_constructible<T>::type* = nullptr>
		TOnce(const char* name = "") : TOnceBase(name), _internal() {}
	
	/* Case 5: Custom collector + name. Delegates to either [1] or [2] */
	template<typename... Ts>
		TOnce(const char* name, void (*fn)(T&, const T&), Ts&&... args) : TOnce(name, std::forward<Ts>(args)...) 
		{ _collector = fn; }
	
	/* Case 5.5: same as previous, but with init list. */
	template<typename U = T, typename std::enable_if<util::has_value_type<U>::value>::type* = nullptr>
		TOnce(const char* name, void (*fn)(T&, const T&), std::initializer_list<typename U::value_type> il) 
		: TOnce(name, il) { _collector = fn; }

	void SetName(std::string name, const char* title = "") override {
		if constexpr(std::is_base_of_v<TNamed, T> || util::has_setname<T>::value) {
			_internal.SetName(name.c_str());
			if constexpr(std::is_base_of_v<TNamed, T>)
				if(title && *title) _internal.SetTitle(title);
		}
        TOnceBase::SetName(name);
	}

	Int_t Write(TFile* f = nullptr, const char* name = "") override {
		if(!name) ERROR("Don't pass nullptr for `name` in TOnce::Write.");
		if(!f || f->IsZombie() || !f->IsOpen()) {
			f = gDirectory->GetFile();
			if(!f || f->IsZombie() || !f->IsOpen())
				ERROR("%s (label: \'%s\') : output ROOT file not supplied or invalid and gDirectory holds no open valid file. (f=0x%016lx)", _name.c_str(), name, (uintptr_t)f);
		}

		if constexpr(std::is_base_of_v<TObject, T>)
			return f->WriteTObject(&_internal, *name ? name : _name.c_str());
		else 
			return f->WriteObject (&_internal, *name ? name : _name.c_str());
	}
    
	void* Load(TFile* f, const char* target = "") override {
		const char* name = *target ? target : _name.c_str();
		if(!name || ! *name) ERROR("Unnamed TOnce<T> object while trying to load from a file.");
		
		if(!f || f->IsZombie() || !f->IsOpen())
			ERROR("(%s) (target: \'%s\') tried to open the file, but didn't receive a handle and also gDirectory holds no open valid file.", _name.c_str(), target);

		auto tmp = std::unique_ptr<T>( f->Get<T>(name) );
		if(!tmp) ERROR("(%s) - asked for " EMPH(%s) " name as key to a static object, got back nullptr.", _name.c_str(), name);
		
		if constexpr(std::is_base_of_v<TH1, T>) {
			tmp->SetDirectory(nullptr);
			_internal.SetDirectory(nullptr);
		}
		
		if constexpr(std::is_base_of_v<TObject, T>) {
			// Try calling `Copy` if it exists in the derived class.
			if constexpr(util::has_copy<T>::value)
				tmp->Copy(_internal);
			else
				_internal = *tmp; // deep-copy ; asserted via type traits on the top.
		}
		else { // Plain STL-types. std::vector, std::array, etc
			_internal = std::move(*tmp);
		}
		
		return (void*)&_internal;
	}
	
	std::unique_ptr<TOnceBase> Clone() const override {
		std::unique_ptr<TOnce<T>> copy;

		if constexpr(std::is_base_of_v<TH1, T>) {
			copy = std::make_unique<TOnce<T>>( this->GetName() );
			this->_internal.Copy( copy->_internal); /* T& -> TH1& -> TObject& */
			copy->_internal.SetDirectory(nullptr);
		}

		else if constexpr(std::is_base_of_v<TObject, T>) {
			copy = std::make_unique<TOnce<T>>( this->GetName() );

			if constexpr(util::has_copy<T>::value)
				this->_internal.Copy(*copy);

			else { /* virtual TObject* Clone(const char *newname="") const */
				if constexpr(util::has_set_directory<T>::value) /* Reassigning the unique_ptr must not double-free. */
					copy->_internal.SetDirectory(nullptr);

				/* This calls case [5] ctor of TOnce<T> */
				copy = std::make_unique<TOnce<T>> (this->GetName(), this->_collector,
						/* T&& */ *static_cast<T*>( _internal.Clone(this->GetName()) )
						);

				if constexpr(util::has_set_directory<T>::value) /* Reassigning the unique_ptr must not double-free. */
					copy->_internal.SetDirectory(nullptr);
			}
		}

		else if constexpr(std::is_copy_constructible_v<T>) {
			copy = std::make_unique<TOnce<T>>(*this);
		}
		else if constexpr(std::is_copy_assignable_v<T>) {
			copy->_internal = _internal;
		}
		else {
			static_assert(! std::is_copy_assignable_v<T>,
				"Wrapped type <T> isn't a ROOT object (inherits from TObject), "
				"and doesn't have copy ctor or copy-assignment operator to clone from!");
		}

		return copy;
	}
 
	void Collect(const TOnceBase& rhs) override {
		const TOnce<T>* cvt = dynamic_cast<const TOnce<T>*> (&rhs);

		if( !cvt) ERROR("Type: %s, wrapped type: %s. TOnce object name \'%s\'. In Collect(..) - dynamic cast failed. "
				"RHS is named \'%s\'", _SELF_TYPE_CSTR, util::type_name<T>().c_str(),
				this->GetName(), rhs.GetName());

		if(strcmp(this->GetName(), rhs.GetName()) != 0)
			ERROR("Type: %s, wrapped type: %s. Trying to collect but objects are called differently. "
					"[1]: %s ; [2]: %s", _SELF_TYPE_CSTR, util::type_name<T>().c_str(),
					this->GetName(), rhs.GetName());

		if(_collector != nullptr)
			return _collector( this->_internal, cvt->_internal );

		/* Go over type traits, try to figure out which call to make, in priority. */
		if constexpr(util::has_dyadic_add_ref<T>::value) {
			return (void)this->_internal.Add( cvt->operator()() );
		}
		else if constexpr(util::has_dyadic_add_ptr<T>::value) {
			return (void)this->_internal.Add( cvt->operator->() );
		}
		else if constexpr(util::has_free_add_fn<T>::value) {
			return (void)Add( this->_internal, cvt->operator()() );
		}
		else if constexpr(util::has_free_mean_fn<T>::value) {
			return (void)Mean( this->_internal, cvt->operator()() );	
		}
		else if constexpr(util::has_add_div_unary_op<T>::value) {
			_internal += cvt->operator()();
			_internal /= 2;
		}
		else {
			ERROR("(%s) - Name: '\%s\' ; Underlying type \'%s\' doesn't define how to add or mean-up two values, "
					"and also its been constructed without a runtime callback. "
					"Define a `void Add(T&, const T& )` function or pass a lambda as second argument "
					"of \'RegisterObject(..)\' (or other ctor).", 
					_SELF_TYPE_CSTR, util::type_name<T>().c_str(), this->GetName());
		}
	}

	T& operator()()             noexcept { return _internal; }
	const T& operator()() const noexcept { return _internal; }

	T* operator->()             noexcept { return &_internal; }
	const T* operator->() const noexcept { return &_internal; }

}; // class TOnce

/* We also here define one free `Add` and `Mean` functions. 
 * Can be specialized/overloaded for certain types. */
template<typename T>
auto Add(T& lhs, const T& rhs) -> decltype((void)(lhs += rhs)) {
	lhs += rhs;
}
template<typename T>
auto Mean(T& lhs, const T& rhs) -> decltype((void)(lhs += rhs, lhs /= 2)) {
	lhs += rhs;
	lhs /= 2;
}

namespace util {
	template<typename T>
	constexpr auto noop_fn() -> void(*)(T&, const T&) {
		return +[](T&, const T&) {};
	}
}
