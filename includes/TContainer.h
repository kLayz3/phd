#pragma once

#include "TObject.h"
#include "TOnce.hxx"
#include <unordered_map>
#include "AuxFunctions.hh"

class TOnceBase;
class TFile;
class TTree;

using TDictInfo = std::unordered_map<std::string, std::string>;
enum class ContainerIO { kINPUT, kOUTPUT }; //!

class TContainer {
	template<typename T> 
	friend class TAnalysisWorker;
	
	std::string _name; //!

protected:
	std::vector<TOnceBase*> _vc; //!
	std::vector<int> _vc_owned;  //!

public:
	TContainer();
	TContainer(std::string );
	virtual ~TContainer();

	virtual void Clean(Option_t* option = "") noexcept;

	inline const char* GetName() const { return _name.c_str(); }
	inline void SetName(std::string(name)) { this->_name = std::move(name); }

	inline void RegisterObjectNoOwn(TOnceBase* b)  { _vc.push_back(b); }
	void RegisterObject(TOnceBase* b);

	template<typename T>
	T* RegisterObject(const char* name, std::initializer_list<typename T::value_type> il) {
		TOnce<T>* obj = new TOnce<T>(name, il);
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		obj->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		RegisterObject(obj);
		
		return obj->operator->(); 
	}

	template<typename T, typename... Ts>
	T* RegisterObject(Ts&&... args) {
		TOnce<T>* obj = new TOnce<T>(std::forward<Ts>(args)...);
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		obj->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		RegisterObject(obj);
		
		return obj->operator->(); 
	}

	/** 
	 * Should be initial (optional) call. Setting up metadata, names and delegating ctors to objects.
	 */
	virtual void Init(TDictInfo info) ;

	/**
	 * Assign the TFile* and TTree* handles to the container.
	 * Switch the container into either input or output mode.
	 * 1. Switch to input mode will load all the static objects and set the branch address of the TTree* to itself, key'ed by its name.
	 * 2. Output mode will create a new TBranch key'ed by its name.
	 */
	virtual void Setup(ContainerIO io_mode, TFile* f = nullptr, TTree* t = nullptr);

	/**
	 * At the end of processing in output mode, writes the output to the TFile* handle.
	 */
	virtual Int_t Write();

	std::vector<TOnceBase*> GetOwnedTOnceObjects() const; //!
	int ClearOwnedTOnceObjects();

	inline TFile* GetFileHandle() const noexcept { return _file_p; } 
	inline TTree* GetTreeHandle() const noexcept { return _tree_p; } 
	
protected:
	ContainerIO _io_mode; //!
	TFile* _file_p = nullptr; //!
	TTree* _tree_p = nullptr; //!

public:
	static TContainer dummy;
	void* _p_self = nullptr; //!

	ClassDef(TContainer, 1);
};

/* This following methods cannot be directly called from base class; every derived class must override it.
 * Problem is that for certain methods, such as:
 * ` template<typename T> TBranch* TTree::Branch(const char* name, T* obj); `
 * This example above above will incorrectly assume type of `this` if it is called from a
 * virtual function from the base class. Even though, `Setup` and `Write` 
 * must be callable via TContainer base class ref/ptr, the derived class must explicitly implement it.
 * That's why all the derived classes must have this boilerplated code... */

#define DECL_CONTAINER_METHODS \
	void Setup(ContainerIO io_mode, TFile* f = nullptr, TTree* t = nullptr) /* override */; \
	Int_t Write() /* override */;


#define IMPL_CONTAINER_METHODS(ClassName) \
	void ClassName::Setup(ContainerIO io_mode, TFile* f, TTree* t) { \
		static_assert(std::is_base_of_v<TContainer, ClassName>); \
		/* Setup the pointers, try to grab defaults if 'f' or 't' are null. */ \
		TContainer::Setup(io_mode, f, t); \
		switch(_io_mode) { \
			case ContainerIO::kINPUT: \
			{ \
				_p_self = (void*)this; \
				Int_t rc = _tree_p ? (_tree_p->SetBranchAddress(this->GetName(), (ClassName**)&_p_self)) : 0; \
				if(rc != 0) \
					ERROR("Container (%s - INPUT) setbranchaddress failed. rc = 0x%08x", GetName(), rc); \
				\
				/* These pointers are set-up either in the ctor or in the `Init` call. */ \
				for(TOnceBase* p : this->_vc) \
					p->Load(_file_p, p->GetName()); \
				break; \
			} \
			case ContainerIO::kOUTPUT: \
			{ \
				if(!_tree_p ) break; \
				/* Container owns the objects that will be written into the ROOT file. */ \
				TBranch* rb = _tree_p->Branch(this->GetName(), this); \
				if(rb == nullptr) \
				ERROR("Container (%s - OUTPUT) creating output branch failed.", GetName()); \
				break; \
			} \
		} \
	} \
	\
	Int_t ClassName::Write() { \
		static_assert(std::is_base_of_v<TContainer, ClassName>); \
		\
		Int_t r = TContainer::Write(); /* Just sanity checks, dead call. */ \
		\
		for(TOnceBase* p : this->_vc) { \
			r += p->Write(_file_p); \
		} \
		/* Writing of the TTree itself is taken care of by the TAnalysisPool instance who ultimately 
		 * hosts all the processes and containers. This class only holds a weak reference to it, to construct
		 * its own branch. */ \
		return r; \
	} \

