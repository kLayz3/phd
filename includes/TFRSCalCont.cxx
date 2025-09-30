#include "TFRSCalCont.h"
#include "AuxFunctions.hh"
#include "TContainer.hxx"
#include "TH2I.h"
#include <string>
#include <tuple>
#include <type_traits>

using nlohmann::json;

json TFRSCalCont::setup{};

void TFRSCalCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", util::type_name<TDictInfo>().c_str());
	auto f = util::get_maybe_ifstream(it->second);
	if(!f.has_value())
		ERROR("File %s not found or not openable.\n", it->second.c_str());
	
	try {
		setup = json::parse(std::move( f.value() ));
	} catch(std::exception const& e) {
		ERROR("Setup parse failed. Reason: %s\n", e.what());
	}

	/* Verify the JSON static + add to param object. */
	for(const auto& [_tpc_i, params] : setup.at("TPC").items()) {
		int i;
		std::string pinfo = params.dump();

		try {
			i = std::stoi(_tpc_i);
		} catch(std::exception const& e) {
			ERROR("TPC \'%s\' unparsable to int. Value: %s\n", _tpc_i.c_str(), pinfo.c_str());
		}
		static const char* keys[] = {
			"x_offset",
			"x_factor",
			"y_offset",
			"y_factor"
		};
		static_assert(std::tuple_size_v<TFRSCalCont::TPCParam> == LEN(keys),
			"Broken parameter mapping construction (array sizes mismatch).\n");

		for(int k=0; k < (int)LEN(keys); ++k) {
			if(! params.at(keys[k]).is_array()) {
				ERROR("TPC%d; Key \'%s\' not found in: %s\n", i, keys[k], pinfo.c_str());	
			}
			try {
				tpc_param->at(k) = params[keys[k]];
			} catch(std::exception const& e) {
				ERROR("Failed setup assignment TPC:%d, keyId:%d, key:%s: reason: %s\n", 
					i, k, keys[k], e.what());
			}
		}
	}
	
	setupName = RegisterObject<std::string>("setup_file", it->second);
	tpc_param = RegisterObject<
		std::remove_reference_t<decltype(*tpc_param)>
	>("tpc_parameters", {});

	h1_ab_s2_before_target = RegisterObject<TH2I>("h1_ab_s2_before_target", "Angle before target (mrad x mrad)", 200, -50, 50, 200, -50, 50);
	h1_xy_s2_before_target = RegisterObject<TH2I>("h1_xy_s2_before_target", "XY before target (mm x mm)", 200, -50, 50, 200, -50, 50);
	h1_xy_s2_after_target = RegisterObject<TH2I>("h1_xy_s2_after_target", "XY after target (mm x mm)", 200, -50, 50, 200, -50, 50);
	h1_ab_s2_after_target = RegisterObject<TH2I>("h1_ab_s2_after_target", "Angle after target (mrad x mrad)", 200, -50, 50, 200, -50, 50);
}

ClassImp(RNFRSCal::Position);
ClassImp(RNFRSCal);
