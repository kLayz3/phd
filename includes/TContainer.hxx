#pragma once

#include "TOnce.hxx"
#include <type_traits>
#include <unordered_map>
#include "AuxFunctions.hh"
#include "ROOT/RNTupleModel.hxx"
#include "TROOT.h"

class TOnceBase;
class TFile;
class TTree;

using TDictInfo = std::unordered_map<std::string, std::string>;
enum class ContainerIO { kINPUT, kOUTPUT }; //!

class TContainerBase {
	inline static std::unique_ptr<ROOT::Experimental::RNTupleModel> _model_read{nullptr};
	inline static std::unique_ptr<ROOT::Experimental::RNTupleModel> _model_write{nullptr};

public:
	/**
	 * Releases the underlying reader model and returns the unique pointer
	 * back to the caller.
	 */
	inline static auto ReleaseModelRead() {
		std::unique_ptr<ROOT::Experimental::RNTupleModel> out = std::move(_model_read);
		// `_model_read` becomes nullptr here.
		return out;
	}
	/**
	 * Releases the underlying writer model and returns the unique pointer
	 * back to the caller.
	 */
	inline static auto ReleaseModelWrite() {
		std::unique_ptr<ROOT::Experimental::RNTupleModel> out = std::move(_model_write);
		// `_model_write` becomes nullptr here.
		return out;
	}
	/**
	 * Create a RNTupleModel only if it's not created already.
	 * Returns the underlying raw pointer.
	 */
	inline static ROOT::Experimental::RNTupleModel* GetModelRead() {
		if(!_model_read) _model_read = ROOT::Experimental::RNTupleModel::Create();
		return _model_read.get();
	}
	/**
	 * Create a RNTupleModel only if it's not created already.
	 * Returns the underlying raw pointer.
	 */
	inline static ROOT::Experimental::RNTupleModel* GetModelWrite() {
		if(!_model_write) _model_write = ROOT::Experimental::RNTupleModel::Create();
		return _model_write.get();
	}
};

template<typename T>
class TContainer : public TContainerBase {
	template<typename U> 
	friend class TAnalysisWorker;

	std::string _filepath;
	std::string _name;
	ContainerIO _io_mode;

	std::vector<TOnceBase*> _vc;
	std::vector<int> _vc_owned;
	
	std::shared_ptr<T> _inner;

public:
	TContainer() {};
	TContainer(std::string name) : _name(std::move(name)) {}
	
	T& operator*() noexcept { return *(_inner.get()); }
	T& inner() noexcept { return this->operator*(); }

	int ClearOwnedTOnceObjects() {
		int r = 0;	
		for(int i : _vc_owned) { delete (_vc.at(i)); ++r; }
		return r;
	}
	virtual ~TContainer() { ClearOwnedTOnceObjects(); }

	void Clean() noexcept {
		if constexpr(util::has_clean_noexcept<T>::value)
			_inner->Clean();
		else
			static_assert(util::has_clean_noexcept<T>::value, "Type has no `void Clean() noexcept` method.\n");
	}

	const char* GetName() const { return _name.c_str(); }
	void SetName(std::string(name)) { this->_name = std::move(name); }

	void RegisterObjectNoOwn(TOnceBase* b)  { _vc.push_back(b); }
	void RegisterObject(TOnceBase* b) {
		_vc_owned.push_back( (int)_vc.size() );
		_vc.push_back(b);
	}

	template<typename U>
	U* RegisterObject(const char* name, std::initializer_list<typename U::value_type> il) {
		TOnce<U>* obj = new TOnce<U>(name, il);
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		obj->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		RegisterObject(obj);
		
		return obj->operator->(); 
	}

	template<typename U, typename... Ts>
	U* RegisterObject(Ts&&... args) {
		TOnce<U>* obj = new TOnce<U>(std::forward<Ts>(args)...);
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		obj->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		RegisterObject(obj);
		
		return obj->operator->(); 
	}

	/** 
	 * Should be initial (optional) call. Setting up metadata, names and delegating ctors to objects.
	 * Override it in the derived class.
	 */
	virtual void Init(TDictInfo info) { (void)info; }

	/**
	 * Assign the TFile* and TTree* handles to the container.
	 * Switch the container into either input or output mode.
	 * 1. Switch to input mode will load all the static objects and set the branch address of the TTree* to itself, key'ed by its name.
	 * 2. Output mode will create a new TBranch key'ed by its name.
	 */
	template<typename Path>
	void Setup(ContainerIO io_mode, Path&& p) {
		static_assert(util::is_pathlike_arg_v<Path>, "Second arg type must be path-like!");
		if(_name.empty())
			ERROR("Unable to proceed with setup call, if the container object is unnamed. "
				EMPH(%s\n), _SELF_TYPE_CSTR);
		if(!util::is_file_readable(std::forward<Path>(p)))
			ERROR("Passed path-like arg (%s) but file instance isn't readable. " EMPH(%s (%s)\n), 
				util::type_name<Path>().c_str(), _name.c_str(), _SELF_TYPE_CSTR);
		auto _fp_maybe = util::get_file_path(std::forward<Path>(p));
		if(!_fp_maybe.has_value())
			ERROR("Unable to fetch file path, did you pass `std::ifstream` by any chance?"
				EMPH(%s (%s)\n), _name.c_str(), _SELF_TYPE_CSTR);

		_filepath = std::move(*_fp_maybe);
		this->_io_mode = io_mode; 

		switch(io_mode) {
			case ContainerIO::kINPUT: 
			{
				auto* m = TContainerBase::GetModelRead();
				_inner = m->MakeField<T>(_name);	
				std::unique_ptr<TFile> f = std::make_unique<TFile>(_filepath.c_str(), "READ");
				for(TOnceBase* p : this->_vc)
					p->Load(f.get(), p->GetName());
			}
			break;
			case ContainerIO::kOUTPUT:
			{
				auto* m = TContainerBase::GetModelWrite();
				_inner = m->MakeField<T>(_name);	
			}
			break;	
		}
	}

	/**
	 * At the end of processing in output mode, write the owned TOnce objects to the ROOT file.
	 */
	Int_t Write() {
		if(_io_mode == ContainerIO::kINPUT)
			ERROR("Bad call to `TContainer::Write()` with its mode being in \'input\'. Name: " 
				EMPH(%s (%s)\n), _name.c_str(), _SELF_TYPE_CSTR);
		if(_filepath.empty())
			ERROR("Call to `Write()` but file path string is empty? Did you call `TContainer::Setup(...)`?"
				EMPH(%s (%s)\n), _name.c_str(), _SELF_TYPE_CSTR);
	
		Int_t r = 0;
		
		std::unique_ptr<TFile> f = std::make_unique<TFile>(_filepath.c_str(), "UPDATE");
		for(TOnceBase* p : this->_vc) {
			r += p->Write(f.get());
		}
		/* Writing of the RNTuple itself is taken care of by the TAnalysisPool instance who ultimately 
		 * hosts all the processes and containers. */ 
		return r;
	}

	std::vector<TOnceBase*> GetOwnedTOnceObjects() const {
		std::vector<TOnceBase*> result{};
		for(auto x : _vc_owned) result.push_back(_vc.at(x));
		return result;
	}
};
