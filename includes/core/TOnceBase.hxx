#pragma once

#include "TH1.h"
class TFile;

struct TOnceBase {
	TOnceBase() { TH1::AddDirectory(kFALSE); }
	TOnceBase(const char* name)  : _name(name) {}
	TOnceBase(std::string name) : _name(std::move(name)) {}

	TOnceBase(const TOnceBase& ) = default;
	TOnceBase& operator=(const TOnceBase& ) = default;
	TOnceBase(TOnceBase&& ) noexcept = default;
	TOnceBase& operator=(TOnceBase&& ) noexcept = default;
	virtual ~TOnceBase() = default; 

	virtual Int_t Write(TFile* file = nullptr, const char* target = "") = 0;
	virtual void* Load(TFile* file, const char* target = "") = 0;
	virtual void Collect(const TOnceBase& ) = 0;
	virtual	std::unique_ptr<TOnceBase> Clone() const = 0;

	inline virtual void SetName(std::string name, const char* title = "") { 
		(void)title;
		_name = std::move(name); 
	}

	inline const char* GetName() const noexcept { return this->_name.c_str(); }

protected:
	std::string _name;
};
