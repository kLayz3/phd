#include "TFRSCalCont.h"
#include "TH2I.h"
#include <stdexcept>
#include <filesystem>

#include "util/JSONParser.h"

using nlohmann::json;
namespace fs = std::filesystem;

TFRSCalCont::TFRSCalCont() : TContainer("FRS") {}

void TFRSCalCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	setup = ParseJSON(file_name);

	/* `file_name` could've been passed relative. Fetch the full path. */
	try {
		fs::path pfull = fs::canonical(file_name);
		setup["file_name"] = pfull.c_str();
	} catch(...) {
		ERROR("Fetching full path from: \'%s\' failed?", file_name.c_str());
	}

	/* Verify the JSON static + add to static param object. */
	for(const auto& [_tpc_i, params] : setup.at("TPC").items()) {
		if(tpc_moniker.find(_tpc_i) == tpc_moniker.end())
			ERROR("TPC parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_tpc_i.c_str(), file_name.c_str(), _tpc_param.size());
		
		u32 i = tpc_moniker.at(_tpc_i);
		if(i >= RNFRSCal::N_VALID_TPC) continue;
		UNROLL_JSON_PARAM(_tpc_param[i], params, 9)
	}

	for(const auto& [_sci_i, params] : setup.at("SCI").items()) {
		if(sci_moniker.find(_sci_i) == sci_moniker.end())
			ERROR("Sci parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_sci_i.c_str(), file_name.c_str(), _sci_param.size());
		
		u32 i = sci_moniker.at(_sci_i);
		if(i >= RNFRSCal::N_VALID_SCI) continue;
		UNROLL_JSON_PARAM(_sci_param[i], params, 3)
	}
}

using T1 = std::remove_reference_t<decltype(*TFRSCalCont::tpc_param)>; // std::array<TPCParam, _>
using T2 = std::remove_reference_t<decltype(*TFRSCalCont::sci_param)>; // std::array<SCIParam, _>

void TFRSCalCont::Setup() {
	std::string setupFileName = "INVALID";
	if(setup.contains("file_name"))
		setupFileName = setup["file_name"].get_ref<const std::string&>();

	for(int i=0; i<RNFRSCal::N_VALID_TPC; ++i) {
		for(int d : {0,1} ) {
			h2_tpc_xy[i][d] = RegisterObject<TH2I>(Form("h2_tpc%d_xy%d", i,d), Form("TPC%d(%d) profile (mm x mm)", i, d), 200, -50, 50, 200, -50, 50);
			h2_tpc_xy[i][d]->GetXaxis()->SetTitle("X [mm]");
			h2_tpc_xy[i][d]->GetYaxis()->SetTitle("Y [mm]");
			h1_tpc_mask[i][d] = RegisterObject<TH1I>(Form("h1_tpc%d_mask%d", i,d), Form("TPC%d(%d) anode mask.", i,d), 4,0,4);

			for(int a: {0,1}) 
				h1_tpc_y[i][2*d+a] = RegisterObject<TH1I>(Form("h1_tpc%d_y%d", i, 2*d+a), Form("TPC%d(%d) anode y-profile (mm)", i, 2*d+a), 800, -40, 40);
		}
	}

	h1_x_sc21_before_target = RegisterObject<TH1I>("h1_x_sc21_before_target", "Position (from Scintillator) before target (mm)", 400,-100,100);
	h1_x_sc22_after_target = RegisterObject<TH1I>("h1_x_sc22_after_target", "Position (from Scintillator) after target (mm)", 400,-100,100);

	/* no-op collector taken deduces from free Add fnc's implemented above. */
	tpc_param = RegisterObject<T1>("tpc_parameters", {});
	sci_param = RegisterObject<T2>("sci_parameters", {});

	/* Just copy over from the static. */
	for(int i=0; i < (int)tpc_param->size(); ++i)
		tpc_param->at(i) = _tpc_param[i];
	for(int i=0; i < (int)sci_param->size(); ++i)
		sci_param->at(i) = _sci_param[i];
	if(SCIParam::channel_to_ns < 0)
		ERROR("Conversion between channel number to ns not given (or parsed) in \'%s\'", setupFileName.c_str());

	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setupFileName);

	h1_x_sc21_before_target->GetXaxis()->SetTitle("X [mm]");
	h1_x_sc22_after_target->GetXaxis()->SetTitle("X [mm]");
}

std::array<double, TPCParam::N_S2_TPC> 
TFRSCalCont::z_s2_tpc(
    std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param
) {
    if(!tpc_param) 
        ERROR("nullptr supplied to z_s2_tpc\n");
    return { 
        tpc_param->at(0).z0,
        tpc_param->at(1).z0,
        tpc_param->at(2).z0,
        tpc_param->at(3).z0
    };
};
std::array<std::array<double, 2>, TPCParam::N_S2_TPC> 
TFRSCalCont::z_s2_tpc_delay_lines(
    std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param
) {
    if(!tpc_param) 
        ERROR("nullptr supplied to z_s2_tpc_delay_lines\n");
    return { 
        tpc_param->at(0).zDL(),
        tpc_param->at(1).zDL(),
        tpc_param->at(2).zDL(),
        tpc_param->at(3).zDL()
    };
}

ClassImp(RNSciCal);
ClassImp(RNSciCal::Measurement);
ClassImp(RNTPCCal);
ClassImp(RNTPCCal::Measurement);
ClassImp(RNFRSCal);

ClassImp(TPCParam);
ClassImp(SCIParam);
