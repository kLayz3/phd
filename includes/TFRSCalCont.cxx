#include "TFRSCalCont.h"
#include "AuxFunctions.hh"
#include "TContainer.hxx"
#include "TH2I.h"
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

using nlohmann::json;

json TFRSCalCont::setup{};

TFRSCalCont::TFRSCalCont() : TContainer("FRS") {
	h1_ab_s2_before_target = RegisterObject<TH2I>("h1_ab_s2_before_target", "Angle before target (mrad x mrad)", 200, -50, 50, 200, -50, 50);
	h1_xy_s2_before_target = RegisterObject<TH2I>("h1_xy_s2_before_target", "XY before target (mm x mm)", 200, -50, 50, 200, -50, 50);
	h1_xy_s2_after_target = RegisterObject<TH2I>("h1_xy_s2_after_target", "XY after target (mm x mm)", 200, -50, 50, 200, -50, 50);
	h1_ab_s2_after_target = RegisterObject<TH2I>("h1_ab_s2_after_target", "Angle after target (mrad x mrad)", 200, -50, 50, 200, -50, 50);
	tpc_param = RegisterObject<
		std::remove_reference_t<decltype(*tpc_param)>
	>("tpc_parameters", {});
}

void TFRSCalCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", util::type_name<TDictInfo>().c_str());
	auto f = util::get_maybe_ifstream(it->second);
	if(!f.has_value())
		ERROR("File \'%s\' not found or not openable.\n", it->second.c_str());
	
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
			if(i < 0 || i >= (int)tpc_param->size())
				throw std::invalid_argument( Form("Arg either negative or parsed as >=%zu", tpc_param->size()) ); 
		} catch(std::exception const& e) {
			ERROR("TPC param: \'%s\' unparsable to int or invalid. Full JSON Value: %s\n"
				"Reason: \'%s\'\n", _tpc_i.c_str(), pinfo.c_str(), e.what());
		}
		static const char* keys[] = {
			"x_offset",
			"x_factor",
			"y_offset",
			"y_factor"
		};
		static_assert(std::tuple_size_v<TFRSCalCont::TPCParam> == LEN(keys),
			"Broken parameter mapping construction (tuple sizes b/w JSON rep and code mismatch).\n");
		
		/* Manually unroll here using a macro. Either that or do aerobatics getting runtime indexing
		 * of tuple (check the pain in TAnalysisPool). */
#define UNROLL_TPC_JSON_PARAM(X) \
		if(! params.at(keys[X]).is_array()) { \
			ERROR("TPC%d; Key \'%s\' not found in: %s\n", i, keys[X], pinfo.c_str());	 \
		} \
		try { \
			std::get<X>( tpc_param->at(i) ) = params[keys[X]]; \
		} catch(std::exception const& e) { \
			ERROR("Failed setup assignment TPC:%d, keyId:%d, key:%s: reason: %s\n",  \
				i, X, keys[X], e.what()); \
		}
		UNROLL_TPC_JSON_PARAM(0)
		UNROLL_TPC_JSON_PARAM(1)
		UNROLL_TPC_JSON_PARAM(2)
		UNROLL_TPC_JSON_PARAM(3)
	}
	
	setupName = RegisterObject<std::string>("setup_file", it->second);
}

ClassImp(RNFRSCal::Position);
ClassImp(RNFRSCal);
