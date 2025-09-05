#include "TContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TROOT.h"
#include <cstdint>
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

void TContainer::Setup(ContainerIO io_mode, TFile* f, TTree* t) {
		/* Maybe users don't want to write the output? */ 
		if(!f || f->IsZombie() || !f->IsOpen()) {
			/* Fine. Try to find it in `gROOT`. */
			std::vector<std::pair<TFile*, TTree*>> _candidate_set{};	
			
			switch(io_mode) {
				case ContainerIO::kINPUT: 
				{
					for(TObject* _f : *gROOT->GetListOfFiles()) {
						f = (TFile*)_f;
						if(f->IsZombie() || !f->IsOpen() || f->IsWritable()) continue;
						for(TObject* _k : *f->GetListOfKeys()) {
							TKey* k = dynamic_cast<TKey*>(_k);
							if(!k) continue;
							TClass* cl = TClass::GetClass(k->GetClassName());
							if(cl && cl->InheritsFrom(TTree::Class())) {
								TTree* _t = dynamic_cast<TTree*>(f->Get(k->GetName()));
								if(!_t || _t->IsZombie()) continue;	
								_candidate_set.emplace_back(f, _t);
							}
						}
					} 
					break;
				}
				/* This is fine since. since if we fetch the object in this call,
				 * and later users (main or TAnalysisPool) want to fetch the object,
				 * they just fetch the same pointer. */
				case ContainerIO::kOUTPUT:
				{
					for(TObject* _f : *gROOT->GetListOfFiles()) {
						f = (TFile*)_f;
						if(f->IsZombie() || !f->IsOpen() || !f->IsWritable()) continue;
						for(TObject* _k : *f->GetList()) {
							TTree* _t = dynamic_cast<TTree*>(_k);
							if(!_t || _t->IsZombie()) continue;
							_candidate_set.emplace_back(f, _t);
						}
					}
					break;	
				}
			}
			if(_candidate_set.size() == 0) 
				ERROR("(%s: \'%s\'), TFile* and TTree* handles not given, and unable to be found inside gROOT.", 
					GetName(), io_mode == ContainerIO::kINPUT ? "Input" : "Output");
			if(_candidate_set.size() > 1) 
				ERROR("(%s: \'%s\'), TFile* and TTree* handles not given, and found multiple TTrees/TFiles readable inside gROOT. Only one deduction is allowed.",
					GetName(), io_mode == ContainerIO::kINPUT ? "Input" : "Output");

			std::tie(f, t) = _candidate_set[0]; // it's asserted, these are not nullptrs.

			/* At this point, we either found correct handles or the call threw. */
			WARN("Deducing (%s: \'%s\'), found TFile*: %p (%s) . TTree*: %p (%s, \'%s\')\n",
				GetName(), io_mode == ContainerIO::kINPUT ? "Input " : "Output", (void*)f, f->GetName(),
				(void*)t, t->GetName(), t->GetTitle());
		}
		/* TFile* passed is non-null and open. 
		 * All file. Maybe users don't want to (de)serialize data row-wise. */
		else if( !t || t->IsZombie()) {
			t = nullptr;
		}

		/* TFile* handle can never be null without exception thrown. */
		this->_file_p  = f; 
		this->_tree_p  = t; 
		this->_io_mode = io_mode; 
}

Int_t TContainer::Write() {
	if(_io_mode == ContainerIO::kINPUT)
		ERROR("Bad call to `TContainer::Write()` with its mode being in \'input\'.");
	if(!_file_p || _file_p->IsZombie() || !_file_p->IsOpen())
		ERROR("TFile* handle which isn't valid or isn't open. Was the TContainer::Setup(...) call successful?");
	return 0;
}

ClassImp(TContainer)
