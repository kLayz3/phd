#pragma once

#include "TOnce.hxx"
#include <unordered_map>
#include "AuxFunctions.hh"
#include "ROOT/RNTupleModel.hxx"
#include "TROOT.h"

class TOnceBase;
class TFile;
class TTree;

using TDictInfo = std::unordered_map<std::string, std::string>;
enum class ContainerIO { kINPUT_FULL, kINPUT_RNONLY, kOUTPUT }; //!

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
		return out;
	}
	/**
	 * Releases the underlying writer model and returns the unique pointer
	 * back to the caller.
	 */
	inline static auto ReleaseModelWrite() {
		std::unique_ptr<ROOT::Experimental::RNTupleModel> out = std::move(_model_write);
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
	
	typedef struct {
		std::string _filepath;
		std::string _name;
		ContainerIO _io_mode;
	} IOInfo;
	/* ^ This is bundled and put behind a single pointer to help L1d locality. 
	 * Usually the data inside is used only during initialisation. */
	std::unique_ptr<IOInfo> _info{};
	std::vector<
		std::unique_ptr<TOnceBase>
	> _vc{};
	
	std::shared_ptr<T> _inner;

public:
	TContainer() : _info( std::make_unique<IOInfo>() ) {}
	TContainer(std::string name) : TContainer() {
		_info->_name = std::move(name);
	}
	
	T& operator*() noexcept { return *(_inner.get()); }
	T& inner() noexcept { return this->operator*(); }

	virtual ~TContainer() = default;

	void Clean() noexcept {
		if constexpr(util::has_clean_noexcept<T>::value)
			_inner->Clean();
		else
			static_assert(util::has_clean_noexcept<T>::value, "Type has no `void Clean() noexcept` method.\n");
	}

	const char* GetName() const { return _info->_name.c_str(); }
	void SetName(std::string(name)) { _info->_name = std::move(name); }

	template<typename U>
	[[nodiscard]] U* RegisterObject(const char* name, std::initializer_list<typename U::value_type> il) {
		std::unique_ptr<TOnce<U>> obj = std::make_unique<TOnce<U>>(name, il);
		TOnce<U>* p = obj.get();
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		p->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		_vc.push_back(std::move(obj));
		
		return p->operator->();
	}

	template<typename U, typename... Ts>
	[[nodiscard]] U* RegisterObject(Ts&&... args) {
		std::unique_ptr<TOnce<U>> obj = std::make_unique<TOnce<U>>(std::forward<Ts>(args)...);
		TOnce<U>* p = obj.get();
		/* Flag the owned object with `CONTAINERNAME_` prefix. */
		p->SetName( ::sstrcat(this->GetName(), "_", obj->GetName()) );
		_vc.push_back(std::move(obj));
		
		return p->operator->();
	}

	/** 
	 * Should be initial (optional) call. Setting up metadata, names and delegating ctors to objects.
	 * Override it in the derived class.
	 */
	virtual void Init(TDictInfo info) { (void)info; }

	/**
	 * Assign the file handles to the container.
	 * Switch the container into either input (2 modes) or output mode.
	 * 1. Switch to `INPUT_FULL` mode will load all the static (TOnce) objects 
	 *    and set the corresponding field of the RNTupleModel to its inner, key'ed by its name.
	 * 2. `INPUT_RNONLY` will only set the field to the RNTupleModel (analog: load the TTree branch).
	 * 3. Output mode `OUTPUT` will create a column entry in the RNTuple key'ed by its name,
	 *    plus keep file name in memory to pass to `Write` in the final TOnce call.
	 */
	template<typename Path>
	void Setup(ContainerIO io_mode, Path&& path) {
		static_assert(util::is_pathlike_arg_v<Path>, "Second arg type must be path-like!");
		if(_info->_name.empty())
			ERROR("Unable to proceed with setup call, if the container object is unnamed. "
				EMPH(%s\n), _SELF_TYPE_CSTR);
	
		/* Check on input modes if the path provided is readable. */
		if(io_mode != ContainerIO::kOUTPUT and !util::is_file_readable(std::forward<Path>(path)))
			ERROR("Passed path-like arg (%s) but file instance isn't readable. " EMPH(%s (%s)\n), 
				util::type_name<Path>().c_str(), _info->_name.c_str(), _SELF_TYPE_CSTR);
		auto _fp_maybe = util::get_file_path(std::forward<Path>(path));
		if(!_fp_maybe.has_value())
			ERROR("Unable to fetch file path, did you pass `std::ifstream` by any chance?"
				EMPH(%s (%s)\n), _info->_name.c_str(), _SELF_TYPE_CSTR);

		_info->_filepath = std::move(*_fp_maybe);
		_info->_io_mode = io_mode; 

		switch(io_mode) {
			case ContainerIO::kINPUT_FULL: {
				std::unique_ptr<TFile> f = std::make_unique<TFile>(
					_info->_filepath.c_str(), "READ");
				for(auto& p : this->_vc)
					p->Load(f.get(), p->GetName());
				[[ fallthrough ]];
			}
			case ContainerIO::kINPUT_RNONLY: {
				auto* m = TContainerBase::GetModelRead();
				_inner = m->MakeField<T>(_info->_name);	
				break;
			}
			case ContainerIO::kOUTPUT: {
				auto* m = TContainerBase::GetModelWrite();
				_inner = m->MakeField<T>(_info->_name);	
				break;	
			}
		}
	}

	/**
	 * At the end of processing in output mode, write the owned TOnce objects to the ROOT file.
	 */
	Int_t Write() {
		if(_info->_io_mode != ContainerIO::kOUTPUT)
			ERROR("Bad call to `TContainer::Write()` with its mode being in \'input\'. Name: " 
				EMPH(%s (%s)\n), _info->_name.c_str(), _SELF_TYPE_CSTR);
		if(_info->_filepath.empty())
			ERROR("Call to `Write()` but file path string is empty? Did you call `TContainer::Setup(...)`?"
				EMPH(%s (%s)\n), _info->_name.c_str(), _SELF_TYPE_CSTR);
	
		Int_t r = 0;
		
		std::unique_ptr<TFile> f = std::make_unique<TFile>(
			_info->_filepath.c_str(), "UPDATE");
		for(const auto& p : _vc) {
			r += p->Write(f.get());
		}
		/* Writing of the RNTuple itself is taken care of by the TAnalysisPool instance who ultimately 
		 * hosts all the processes and containers. */ 
		return r;
	}

	inline std::size_t GetNOwnedTOnceObjects() { return _vc.size(); }

};
