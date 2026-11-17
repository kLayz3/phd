#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TPad.h>
#include <TROOT.h>
#include <cstring>
#include <memory>
#include <string>
#include <variant>

void extract_pad(
	int padn,
	std::variant<
		TFile*, std::string
	> source = gFile,
	const char* cname = nullptr
) {
	TDirectory::TContext directory_context;
	std::unique_ptr<TFile> owned_file;
	TFile* f = nullptr;
	if(auto* pointer = std::get_if<TFile*>(&source)) {
		f = *pointer;
	} else {
		const auto& filename = std::get<std::string>(source);
		owned_file.reset(TFile::Open(filename.c_str(), "READ"));
		f = owned_file.get();
	}
	if(!f || f->IsZombie()) {
		Error("extract_pad", "No readable ROOT file supplied or currently open");
		return;
	}
	
	TCanvas* c = nullptr;
	if(!cname || !strlen(cname)) {
		// Try to find a canvas object inside ROOT.
		for(auto* object: *f->GetListOfKeys()) {
			const auto* key = static_cast<TKey*>(object);
			auto* type = TClass::GetClass(key->GetClassName());
			if(!type || !type->InheritsFrom(TCanvas::Class())) continue;
			c = dynamic_cast<TCanvas*>(f->Get(key->GetName()));
			if(c) break;
		}
	}
	else {
		c = dynamic_cast<TCanvas*>(f->Get(cname));
	}
	
	if(!c) return;

	auto* pad = dynamic_cast<TPad*>(c->GetPad(padn));
	if(!pad) return;

	auto* out = new TCanvas("out", "Extracted pad", 1600, 1000);
	out->cd();

	// Cloned histograms must outlive a file opened just for this invocation.
	TDirectory::TContext clone_context(gROOT);
	auto* cpy = static_cast<TPad*>(pad->Clone("padc"));
	cpy->Draw();
	cpy->SetPad(0.,0., 1.,1.);
	cpy->Modified();

	out->Modified();
	out->Update();
}
