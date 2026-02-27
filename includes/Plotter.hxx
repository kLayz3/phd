/* A small helper library to wrap around different pretty-plot styles for histograms, etc. 
 * These are called in macros, and don't expect to be rather optimized.
 */

#pragma once
#include <string>
#include <cstring>
#include "TH1.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TNamed.h"

namespace mnd {

inline std::string fname_short(std::string fileName) {
	auto fname_short = std::string(strrchr(fileName.c_str(), '/') ? strrchr(fileName.c_str(), '/') + 1 : fileName); 
	return fname_short.substr(0, fname_short.find(".root"));
}

inline void SetTitle(TObject* x, std::string title) {
	if(!x) throw std::runtime_error("SetTitle: null TObject*");
	if(auto* n = dynamic_cast<TNamed*>(x)) {
		n->SetTitle(title.c_str()); return;
	}
	if(auto* c = dynamic_cast<TCanvas*>(x)) {
		c->SetTitle(title.c_str()); return;
	}
	throw std::runtime_error(
		std::string("Type '") + x->ClassName() + "' doesn't have a supported SetTitle overload"
	);
}

}
