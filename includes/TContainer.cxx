#include "TContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include <cstring>

TContainer::TContainer() {}
TContainer::TContainer(std::string name) : _name(std::move(name)) {}
TContainer::~TContainer() { ClearOwnedTOnceObjects(); }

void TContainer::Init(TDictInfo info) { (void)info; }
void TContainer::Clean(Option_t* option) noexcept { (void)option; }

TContainer TContainer::dummy{};

void TContainer::RegisterObject(TOnceBase* b) {
	_vc_owned.push_back( (int)_vc.size() );
	_vc.push_back(b);
}

std::vector<TOnceBase*> TContainer::GetOwnedTOnceObjects() const {
	std::vector<TOnceBase*> result{};
	for(auto x : _vc_owned) result.push_back(_vc.at(x));
	return result;
}

int TContainer::ClearOwnedTOnceObjects() {
	int r = 0;	
	for(int i : _vc_owned) { delete (_vc.at(i)); ++r; }
	return r;
}

void TContainer::Setup(TFile* f, TTree* t, ContainerIO io_mode) {
		/* Maybe users don't want to write the output? */ 
		if(!f || f->IsZombie() || !f->IsOpen()) { 
			WARN("Container (%s - %s) - passed TFile* handle which: isn't in gDirectory or pointer isn't valid, isn't opened, could be ignored. Is OK.", GetName(), (io_mode == ContainerIO::kINPUT) ? "INPUT" : "OUTPUT" ); 
			return; 
		} 
 
		if(!t || t->IsZombie()) { 
			WARN("Container (%s - %s) - passed bad TTree handle to the function. Ignoring it, is OK. Event-by-event data won't be written/read.", GetName(), (io_mode == ContainerIO::kINPUT) ? "INPUT" : "OUTPUT"); 
		} else { 
			this->_tree_p = t; 
		} 
 
		this->_io_mode = io_mode; 
		this->_file_p = f; 
}

Int_t TContainer::Write() {
	if(!_file_p || _file_p->IsZombie() || !_file_p->IsOpen())
		ERROR("TFile* handle which isn't valid or isn't open. Was the TContainer::Setup call successful?");
	if(!_tree_p || _tree_p->IsZombie())
		ERROR("TTree* handle which isn't valid or isn't open. Was the TContainer::Setup call successful?");
	if(_io_mode == ContainerIO::kINPUT)
		ERROR("Forbidden to call TContainer::Write with mode switch to input instead of output/mixed.");
	return 0;
}


ClassImp(TContainer)
