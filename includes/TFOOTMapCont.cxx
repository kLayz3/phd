#include "TFOOTMapCont.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TGraph.h"
#include <cmath>
#include <string>

TFOOTMapCont::TFOOTMapCont(int N) : TContainer(Form("FOOT%d", N)), FOOT_N(N) {}

void TFOOTMapCont::Init(TDictInfo info) {
	auto n_it = info.find("FOOT_ID");
	if(n_it == info.end())
		ERROR("FOOT_ID key not found in the info hashmap.");
	i32 n;
	try {
		n = std::stoi(n_it->second);
	} catch(const std::exception& e) {
		ERROR("FOOT_ID found: " EMPH(%s) " but unparsable to integer. Err: %s", n_it->second.c_str(), e.what());  
	}
	FOOT_N = n;
	this->SetName(Form("FOOT%d", n));
}

ClassImp(RNFOOTMap);
