#include "TContainer.h"
#include "TOnce.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include <cstring>
#include "AuxFunctions.hh"

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

template<typename T>
T* TContainer::RegisterObject(const char* name, std::initializer_list<typename T::value_type> il) {
	TOnce<T>* obj = new TOnce<T>(name, il);
	/* Flag the owned object with `CONTAINERNAME_` prefix. */
	obj->SetName( sstrcat(this->GetName(), "_", obj->GetName()) );
	RegisterObject(obj);
	
	return obj->operator->(); 
}

template<typename T, typename... Ts>
T* TContainer::RegisterObject(Ts&&... args) {
	TOnce<T>* obj = new TOnce<T>(std::forward<Ts>(args)...);
	/* Flag the owned object with `CONTAINERNAME_` prefix. */
	obj->SetName( sstrcat(this->GetName(), "_", obj->GetName()) );
	RegisterObject(obj);
	
	return obj->operator->(); 
}

void TContainer::Write(TFile* f, TTree* t) {
	if(!f) f = gDirectory->GetFile();	
	if(!f || f->IsZombie() || !f->IsOpen())
		ERROR("Passed TFile* handle which isn't valid or isn't open.");
	
	if(!t) { // Try find it in the upper most directory.
		TIter next(f->GetListOfKeys());
		TKey* key;
		while((key = (TKey*)next()) != nullptr) {
			t = dynamic_cast<TTree*>( key->ReadObj() );
			if(t && !t->IsZombie())
				return this->Write(f,t); // take first TTree* you find.
		}
	}
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

IMPL_CONTAINER_SETUP(TContainer);
ClassImp(TContainer)
