#pragma once
#include "Rtypes.h"
#include "TH1.h"

class TFile;

class TOnceBase {
protected:
	std::string _name;

public:
	TOnceBase() { TH1::AddDirectory(kFALSE); }
	TOnceBase(const char* name) : TOnceBase() { _name = name; }
	TOnceBase(std::string name) : TOnceBase() { _name = std::move(name); }

	TOnceBase(const TOnceBase& ) = default;
	TOnceBase& operator=(const TOnceBase& ) = default;
	TOnceBase(TOnceBase&& ) noexcept = default;
	TOnceBase& operator=(TOnceBase&& ) noexcept = default;
	virtual ~TOnceBase() = default; 

	virtual Int_t Write(TFile* file = nullptr, const char* target = "") = 0;
	virtual void* Load(TFile* file, const char* target = "") = 0;

	virtual void SetName(std::string name, const char* title = "");
	const char* GetName() const;
};
