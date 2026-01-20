#include "TFRSCalCont.h"
#include "TH2I.h"
#include <stdexcept>

using nlohmann::json;

TFRSCalCont::TFRSCalCont() : TContainer("FRS") {}

void TFRSCalCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	auto f = mnd::get_maybe_ifstream(file_name);
	if(!f.has_value())
		ERROR("File \'%s\' not found or not openable.\n", file_name.c_str());
	
	try {
		setup = json::parse(std::move( f.value() ));
	} catch(std::exception const& e) {
		ERROR("Setup parse failed. Reason: %s\n", e.what());
	}
	setup["file_name"] = file_name;

	/* Verify the JSON static + add to static param object. */
	for(const auto& [_tpc_i, params] : setup.at("TPC").items()) {

		if(tpc_moniker.find(_tpc_i) == tpc_moniker.end())
			ERROR("TPC parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_tpc_i.c_str(), file_name.c_str(), _tpc_param.size());
		
		std::string pinfo = params.dump();
		u32 i = tpc_moniker.at(_tpc_i);
		if(i >= RNFRSCal::N_VALID_TPC) continue;

		static const char* keys[] = {
			"x_offset",
			"x_factor",
			"y_offset",
			"y_factor",
			"csum_lim",
			"anode_diff_lim",
			"dl_left_diff_lim",
			"dl_right_diff_lim",
			"sci_ref_lim",
			"z0"
		};
		static_assert(TPCParam::N_PARAMS == mnd::len(keys),
			"Broken parameter mapping construction (tuple sizes b/w JSON representation and code mismatch)\n");
		
		/* Manually unroll here using a macro. Either that or do aerobatics getting runtime indexing. */
#define UNROLL_TPC_JSON_PARAM(X) \
		try { \
			/* Hacky part: right side is explicit conversion from nlohmann::json object
			 * to std::array<?,?>. Type traits explore the exact nested type. */ \
			_tpc_param.at(i).get<X>() = params[keys[X]] \
				.get< \
					std::remove_reference_t< \
						decltype( std::declval<TPCParam&>().get<X>() ) \
					> \
				>(); \
		} catch(std::exception const& e) { \
			ERROR("Failed setup assignment TPC:%d, keyId:%d, key:%s: reason: %s\n",  \
				i, X, keys[X], e.what()); \
		}

		UNROLL_TPC_JSON_PARAM(0)
		UNROLL_TPC_JSON_PARAM(1)
		UNROLL_TPC_JSON_PARAM(2)
		UNROLL_TPC_JSON_PARAM(3)
		UNROLL_TPC_JSON_PARAM(4)
		UNROLL_TPC_JSON_PARAM(5)
		UNROLL_TPC_JSON_PARAM(6)
		UNROLL_TPC_JSON_PARAM(7)
		UNROLL_TPC_JSON_PARAM(8)
		UNROLL_TPC_JSON_PARAM(9)
	}

	for(const auto& [_sci_i, params] : setup.at("SCI").items()) {
		if(sci_moniker.find(_sci_i) == sci_moniker.end())
			ERROR("Sci parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_sci_i.c_str(), file_name.c_str(), _sci_param.size());
		
		std::string pinfo = params.dump();
		u32 i = sci_moniker.at(_sci_i);
		if(i >= RNFRSCal::N_VALID_SCI) continue;
		
		static const char* keys[] = {
			"x_offset",
			"x_factor",
			"cdiff_lim",
			"z0"
		};
		static_assert(SCIParam::N_PARAMS == mnd::len(keys),
			"Broken parameter mapping construction (tuple sizes b/w JSON representation and code mismatch)\n");

		/* Copy-pasta from above. */
#define UNROLL_SCI_JSON_PARAM(X) \
		try { \
			_sci_param.at(i).get<X>() = params[keys[X]] \
				.get< \
					std::remove_reference_t< \
						decltype( std::declval<SCIParam&>().get<X>() ) \
					> \
				>(); \
		} catch(std::exception const& e) { \
			ERROR("Failed setup assignment SCI%s (%u), key: \'%s\', reason: %s\n", \
				_sci_i.c_str(), i, #X, e.what()); \
		}
		
		UNROLL_SCI_JSON_PARAM(0)
		UNROLL_SCI_JSON_PARAM(1)
		UNROLL_SCI_JSON_PARAM(2)
		UNROLL_SCI_JSON_PARAM(3)
	}
}

/* A small Add function for this type, which is a no-op. If it's not defined,
 * then folding `Collect` over this type will be compile error. */
using T1 = std::remove_reference_t<decltype(TFRSCalCont::_tpc_param)>; // std::array<TPCParam, _>
using T2 = std::remove_reference_t<decltype(TFRSCalCont::_sci_param)>; // std::array<SCIParam, _>
void Add(T1&, const T1&) {}
void Add(T2&, const T2&) {}

void TFRSCalCont::Setup() {
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
		ERROR("Conversion between channel number to ns not given (or parsed) in \'%s\'", setup["file_name"].get_ref<const std::string&>().c_str());

	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());

	h1_x_sc21_before_target->GetXaxis()->SetTitle("X [mm]");
	h1_x_sc22_after_target->GetXaxis()->SetTitle("X [mm]");
}

ClassImp(RNSciCal);
ClassImp(RNSciCal::Measurement);
ClassImp(RNTPCCal);
ClassImp(RNTPCCal::Measurement);
ClassImp(RNFRSCal);

ClassImp(TPCParam);
ClassImp(SCIParam);
