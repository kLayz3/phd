#pragma once

#include "TObject.h"
#include <unordered_map>

class TOnceBase;
class TFile;
class TTree;
class TAnalysisWorker;

using TDictInfo = std::unordered_map<std::string, std::string>;
enum class ContainerIO { kINPUT, kOUTPUT }; //!

class TContainer {
	friend class TAnalysisWorker; //!
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
	T* RegisterObject(const char*, std::initializer_list<typename T::value_type>);
	template<typename T, typename... Ts>
	T* RegisterObject(Ts&&... args);

	/** 
	 * Should be initial (optional) call. Setting up the names and objects.
	 */
	virtual void Init(TDictInfo info) ;

	/**
	 * Call after that to switch the container into either input or output mode.
	 */
	virtual void Setup(TFile* f = nullptr, TTree* t = nullptr, ContainerIO io_mode = ContainerIO::kINPUT);

	/**
	 * At the end of processing in output mode, this call writes the output to
	 * the root file.
	 */
	void Write(TFile* f = nullptr, TTree* t = nullptr);

	std::vector<TOnceBase*> GetOwnedTOnceObjects() const;
	int ClearOwnedTOnceObjects();

protected:
	ContainerIO _io_mode; //!
	TFile* _file_p = nullptr; //!
	TTree* _tree_p = nullptr; //!

public:
	static TContainer dummy;

	ClassDef(TContainer, 1);
};

/* This following method cannot be kept as a dummy, every derived class must override it.
 * Problem is that for certain methods, such as:
 * ` template<typename T> TBranch* TTree::Branch(const char* name, T* obj); `
 * This function above will incorrectly assume type of `this` if it is called from a non-overridden 
 * virtual function from the base class. Though, `Setup` must be callable via TContainer base class ref/ptr.
 * That's why all the derived classes must override it with the boilerplated code... */

#define IMPL_CONTAINER_SETUP(ClassName) \
	void ClassName::Setup(TFile* f, TTree* t, ContainerIO io_mode) { \
		static_assert(std::is_base_of_v<TContainer, ClassName>); \
		/* Maybe users don't want to write the output? */ \
		if(!f || f->IsZombie() || !f->IsOpen()) { \
			WARN("Container (%s - %s) - passed TFile* handle which: isn't in gDirectory or pointer isn't valid, isn't opened, could be ignored. Is OK.", GetName(), (io_mode == ContainerIO::kINPUT) ? "INPUT" : "OUTPUT" ); \
			return; \
		} \
 \
		if(!t || t->IsZombie()) { \
			WARN("Container (%s - %s) - passed bad TTree handle to the function. Ignoring it, is OK. Event-by-event data won't be written/read.", GetName(), (io_mode == ContainerIO::kINPUT) ? "INPUT" : "OUTPUT"); \
		} else { \
			this->_tree_p = t; \
		} \
 \
		this->_io_mode = io_mode; \
		this->_file_p = f; \
		\
		/* Set the instance to be either the input or the output. */  \
		switch(io_mode) { \
			case ContainerIO::kINPUT: { \
					Int_t rc = t->SetBranchAddress(this->GetName(), this); \
					if(rc != 0) \
					ERROR("Container (%s - INPUT) setbranchaddress failed. rc = 0x%08x", GetName(), rc);	 \
					\
					/* These pointers are set-up either in the ctor or in the `Init` call. */ \
					for(TOnceBase* p : this->_vc) \
					p->Load(f, sstrcat(this->GetName(), "_", p->GetName()).c_str()); \
				} \
				break; \
			case ContainerIO::kOUTPUT: { \
					if(!t || t->IsZombie()) break; \
					/* Container owns the objects that will be written into the ROOTFILE */ \
					TBranch* rb = t->Branch(this->GetName(), this); \
					if(rb == nullptr)  \
					ERROR("Container (%s - OUTPUT) creating output branch failed.", GetName()); \
				} \
				break; \
		} \
	} \

