#pragma once

#include <filesystem>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "TROOT.h"
#include "TCanvas.h"
#include "TInterpreter.h"

using A2 = std::array<double, 2>;
using A3 = std::array<double, 3>;

enum class DoSave { yes, no };

namespace canvas {
	enum struct Extension { png, jpeg, pdf, C, root };

	inline void save_all (
		Extension extension = Extension::png, 
		std::vector<std::string_view> extra_tag = {}
	) { 
		const char* ext;
#define HANDLE_CASE_SAVE_ALL(e) case(Extension::e): { ext = #e; break; }

		switch(extension) {
			HANDLE_CASE_SAVE_ALL(png)
			HANDLE_CASE_SAVE_ALL(jpeg)
			HANDLE_CASE_SAVE_ALL(pdf)
			HANDLE_CASE_SAVE_ALL(C)
			HANDLE_CASE_SAVE_ALL(root)
		}
#undef HANDLE_CASE_SAVE_ALL

		const char* macro_name = gInterpreter->GetCurrentMacroName(); 
		std::string stem = std::filesystem::path(macro_name).stem().string();

		std::filesystem::path p = "autosave";
		p = p / stem / ext;
		for(auto& tag: extra_tag)
			p /= tag;

		try {
			std::filesystem::create_directories(p); 
		} catch(const std::filesystem::filesystem_error& e) {
			fprintf(stderr, "canvas::save_all : Error in creating directories: \'%s\', err: %s\n", p.c_str(), e.what());  
			return;
		}

		std::vector<TCanvas*> cs;
		for(TObject* k_ : *gROOT->GetListOfCanvases()) {
			if(TCanvas* c = dynamic_cast<TCanvas*>(k_))
				cs.push_back(c);
		}
		if(cs.empty()) return;

		for(auto* c : cs) {
			auto name = std::string(c->GetName()); 
			auto outfile = p / (name + "." + ext);
			// Force rendering
			c->Modified();
			c->Update();
			c->SaveAs( outfile.c_str() ); 
		}
		if(cs.size() == 1 or (extension != Extension::png and extension != Extension::jpeg)) return;

		// In case of two or more canvases saved, also collect them into a .pdf
		std::string outpdf = (p / "all.pdf").string();

		cs.front()->Print(Form("%s(", outpdf.c_str()));
		for(size_t i=1; i < cs.size() - 1; ++i) {
			cs[i]->Print( outpdf.c_str() );
		} 
		cs.back()->Print(Form("%s)", outpdf.c_str()));
	}
};
