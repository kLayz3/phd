#pragma once

#include "AuxFunctions.hh"
#include <unordered_map>
#include "RtypesCore.h"
#include "TOnce.hxx"

#include "TBuffer.h"
#include "TBufferFile.h"

using TDictInfo = std::unordered_map<std::string, std::string>;

struct TContainerBase {
	std::string _name;
	
	inline const char* GetName() const noexcept { return _name.c_str(); }
	inline void SetName(std::string(name)) noexcept { _name = std::move(name); }

	TContainerBase() = default;
	TContainerBase(std::string name) : _name(std::move(name)) {}

	/** 
	 * Should be initial (optional) call. Setting up persistent metadata, and names.
	 */
	virtual void Init(TDictInfo info) { (void)info; }

	/**
	 * Delayed construction of TOnce objects. Possibly override it in the derived classes.
	 */
	virtual void Setup() = 0;
};

struct MTDeepCopy { explicit MTDeepCopy() = default; };

/**
 * Encapsulates a type needed to be used as a (de)serialization target
 * from/to RNTuple column. It still is an abstract type, since the `Setup()` method isn't defined.
 * Users are at most expected to define the Setup() method, where the `RegisterObject` methods will be
 * chained, to give back the raw resource handles to the users' main container class.
 */
template<typename T>
struct TContainer : TContainerBase {
	//static_assert( std::is_copy_constructible_v<T>, "Type must be copy constructible.");
	//static_assert( std::is_copy_assignable_v<T>, "Type must be copy assignable.");
	//static_assert( std::is_move_constructible_v<T>, "Type must be move constructible.");
	//static_assert( std::is_move_assignable_v<T>, "Type must be move assignable.");
	static_assert(!std::is_void_v<T>, "Hello?");
	static_assert(!std::is_array_v<T>, "Must not pass raw C-style arrays. Pass an `std::array<T,N>` instead.");
	static_assert(!std::is_pointer_v<T>, "Must not pass pointer type (T*).");
	static_assert(!std::is_reference_v<T>, "Must not pass ref type (T&).");

	template<typename>    friend struct TProcessor;
	template<typename...> friend struct TAnalysisProcess;

	using inner_type = T;
	using TOnceBaseVec = std::vector <
		std::shared_ptr<TOnceBase>
	>;

protected:
	TOnceBaseVec _vc;

private:
	std::shared_ptr<T> _inner;

public:
	TContainer() {};
	TContainer(std::string name) : TContainerBase(name) {}

	/* Ok, some clarification. Lets say `struct Derived : TContainer<T>` is the child class.
	 * Then the next two methods below must get called in `Derived::Setup()`
	 * We create the objects, pulling the read objects from disk, and receiving back a shared handle of the resource. 
	 * This works fine for the original TContainer<T> instance.
	 * But when cloning a TAnalysisProcess, it clones the underlying 'TProcessor<Out(Ins...)>' types. 
	 * Instance behind *this* pointer, for the reader, is created via the copy ctor of the TContainer<T>,
	 * while the write instances are default constructed. The clone's ctor will also call this sequence initially, each pulling its own
	 * copy of the read-objects from disk into RAM. Using unnecessary disk space, if simply all the clones own a 
	 * read-only shared pointer to a single instance of such objects.
	 * - key difference is that then the Original singleton will simply switch the clones' variant to non-owning
	 * type, and the temporarily created (unique) objects of the clone will get deleted. 
	 * This is important for read containers - as only one set of these `TOnce<T>` objects will get (de)serialized.
	 * Clones will only hold a raw pointer vector variant, and the original Container holds the owning pointers. 
	 * For write containers, each clone has its own unique copy. */

	/**
	 * Create an object to be serialized to/from a ROOT file. Two overloads exist for passing an
	 * initializer list ctor (like for `std::array<T>`) or a basic forwarding ctor (emplace-style). 
	 */
	template<typename U>
	[[nodiscard]] U* RegisterObject(const char* name, std::initializer_list<typename U::value_type> il) {
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		std::string obj_name = util::sstrcat(this->GetName(), "_", name); 
		for(auto& o : _vc) {
			if(! strcmp(o->GetName(), obj_name.c_str())) {
				TOnce<U> *dcast = dynamic_cast<TOnce<U>*>( o.get() );
				if(!dcast)
					ERROR("Attempted to retrieve object named \'%s\' from %zu-sized list of owned TOnce<..> objects."
						"Found at address 0x%lx a 'TOnceBase' but dynamic_cast failed? Same object registered multiple times? (%s)",
						name, _vc.size(), (uintptr_t)o.get(), _SELF_TYPE_CSTR);
				return dcast->operator->();
			};
		}
		
		std::shared_ptr<TOnce<U>> obj = std::make_shared<TOnce<U>>(obj_name.c_str(), il);
		TOnce<U>* p = obj.get();

		_vc.push_back( std::move(obj) );
		return p->operator->();
	}

	template<typename U, typename... Ts>
	[[nodiscard]] U* RegisterObject(const char* name, Ts&&... args) {
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		std::string obj_name = util::sstrcat(this->GetName(), "_", name); 
		for(auto& o : _vc) {
			if(! strcmp(o->GetName(), obj_name.c_str())) {
				TOnce<U> *dcast = dynamic_cast<TOnce<U>*>( o.get() );
				if(!dcast)
					ERROR("Attempted to retrieve object named \'%s\' from %zu-sized list of owned TOnce<..> objects."
						"Found at address 0x%lx a 'TOnceBase' but dynamic_cast failed? Same object registered multiple times? (%s)",
						name, _vc.size(), (uintptr_t)o.get(), _SELF_TYPE_CSTR);
				return dcast->operator->();
			};
		}
		
		std::shared_ptr<TOnce<U>> obj = std::make_shared<TOnce<U>>(obj_name.c_str(), std::forward<Ts>(args)...);
		TOnce<U>* p = obj.get();

		_vc.push_back( std::move(obj) );
		return p->operator->();
	}

	void Clean() noexcept {
		if constexpr(util::has_clean_noexcept<T>::value)
			_inner->Clean();
		else
			static_assert(util::has_clean_noexcept<T>::value, "Type <T> has no `void Clean() noexcept` method.\n");
	}

	/* Rule of five. We keep it here to be explicit about copy-ctor. */
	TContainer(const TContainer& )                = default;
	TContainer& operator=(const TContainer& )     = default;
	TContainer(TContainer&& )            noexcept = default;
	TContainer& operator=(TContainer&& ) noexcept = default;
	~TContainer()                                 = default;
	
	const TOnceBaseVec& GetTOnceVec() noexcept { return this->_vc; }

	T& operator*()             noexcept { return *_inner; }
	const T& operator*() const noexcept { return *_inner; }
	T& inner()                 noexcept	{ return this->operator*(); }
	const T& inner()     const noexcept	{ return this->operator*(); }
	T* raw()                   noexcept	{ return _inner.get(); }
	const T* raw()       const noexcept	{ return _inner.get(); }
};

/* Type `T` is either something like TXXXYYYEvent (Go4) or a custom structure
 * for UCESB. 
 * Namely, UCESB splits all the trees directly into leaves, so address has to be mapped
 * sequentally.
 * Encapsulates a type needed to be used as a deserialization target from a TTree branch,
 * that ended up as raw-baked into an RNTuple. Once an object gets materialized and we get
 * back a byte-array, uset the Streamer to materialize the concrete object.
 */
template<typename T>
struct TRawContainer : TContainerBase {
	static_assert(!std::is_void_v<T>, "Hello?");
	static_assert(!util::is_an_array_v<T>, "Must not pass raw C-style array or `std::array<T,N>`.");
	static_assert(!std::is_pointer_v<T>, "Must not pass pointer type (T*).");
	static_assert(!std::is_reference_v<T>, "Must not pass ref type (T&).");
	
	template<typename>    friend struct TProcessor;
	template<typename...> friend struct TAnalysisProcess;
	
	//using TRaw = std::vector<unsigned char>;
	//using inner_type = TRaw;
	using inner_type = T;

private:
	//std::shared_ptr <
	//	std::vector<unsigned char>
	//> _inner;
	//std::shared_ptr<T> _inner_raw;
	T* _inner;

public:
	TRawContainer() = default;
	TRawContainer(std::string name) : TContainerBase(name) {}

	void Setup() override {}

	//void Materialize() {
	//	TBufferFile buf(TBuffer::kRead, _inner->size(),
	//		_inner->data(), /*own=*/kFALSE);

	//	T* obj = static_cast<T*>( buf.ReadObject(T::Class()) );
	//	_inner.reset(obj);
	//}

	T& operator*()             noexcept { return *_inner; }
	const T& operator*() const noexcept { return *_inner; }
	T& inner()                 noexcept	{ return this->operator*(); }
	const T& inner()     const noexcept	{ return this->operator*(); }
	T* raw()                   noexcept	{ return _inner; }
	const T* raw()       const noexcept	{ return _inner; }
};
