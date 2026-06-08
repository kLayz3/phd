#include "TFRSHitCont.h"
#include "util/JSONParser.h"

using nlohmann::json;

static nlohmann::json setup {};
static FRSIdParam _s2p, _s3p;
static FRSTargetParam _sTar;
static std::array<TPCParam, RNFRSCal::N_VALID_TPC> _tpc_param {};
static std::array<SCIParam, RNFRSCal::N_VALID_SCI> _sci_param {};

TFRSHitCont::TFRSHitCont() : TContainer("FRS") {}  
 
void TFRSHitCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;

	setup = ParseJSON(file_name);
	setup["file_name"] = file_name;

	/* Copy over params from cal step. */
	for(const auto& [_tpc_i, params] : setup.at("TPC").items()) {
		if(TFRSCalCont::tpc_moniker.find(_tpc_i) == TFRSCalCont::tpc_moniker.end())
			ERROR("TPC parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_tpc_i.c_str(), file_name.c_str(), _tpc_param.size());
		
		u32 i = TFRSCalCont::tpc_moniker.at(_tpc_i);
		if(i >= RNFRSCal::N_VALID_TPC) continue;
		UNROLL_JSON_PARAM(_tpc_param[i], params, 9)
	}
	for(const auto& [_sci_i, params] : setup.at("SCI").items()) {
		if(TFRSCalCont::sci_moniker.find(_sci_i) == TFRSCalCont::sci_moniker.end())
			ERROR("Sci parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_sci_i.c_str(), file_name.c_str(), _sci_param.size());
		
		u32 i = TFRSCalCont::sci_moniker.at(_sci_i);
		if(i >= RNFRSCal::N_VALID_SCI) continue;
		UNROLL_JSON_PARAM(_sci_param[i], params, 3)
	}

	auto j_it = setup.find("FRS");
	if(j_it == setup.end()) ERROR("\'FRS\' key not found in JSON file: %s\n", file_name.c_str());
	auto& jfrs = j_it.value();

	j_it  = jfrs.find("S2");
	if(j_it == jfrs.end()) ERROR("\'S2\' key not found in \'FRS\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_s2p, j_it.value(), 5);
	
	/* Unpack the S2 target param values. */
	auto& js2 = j_it.value();
	j_it = j_it.value().find("target");
	if(j_it == js2.end()) ERROR("\'targte\' key not found in \'S2\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_sTar, j_it.value(), 1);
	
	j_it  = jfrs.find("S3");
	if(j_it == jfrs.end()) ERROR("\'S3\' key not found in \'FRS\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_s3p, j_it.value(), 5);
}

void Add(FRSIdParam&, const FRSIdParam&) {}
void Add(FRSTargetParam&, const FRSTargetParam&) {}

void TFRSHitCont::Setup() {
	h2_track_x = RegisterObject<TH2D>("h2_s2_track_x", "S2 Tracking;z[mm];x[mm]", 400, 0, 4200, 200, -100, 100);
	h2_track_y = RegisterObject<TH2D>("h2_s2_track_y", "S2 Tracking;z[mm];y[mm]", 400, 0, 4200, 200, -100, 100);
	/* no-op collector taken deduces from free Add fnc's implemented above. */
	tpc_param = RegisterObject<std::array<TPCParam, RNFRSCal::N_VALID_TPC>>("tpc_parameters", {});
	sci_param = RegisterObject<std::array<SCIParam, RNFRSCal::N_VALID_SCI>>("sci_parameters", {});

	/* Just copy over from the static. */
	for(int i=0; i < (int)tpc_param->size(); ++i)
		tpc_param->at(i) = _tpc_param[i];
	for(int i=0; i < (int)sci_param->size(); ++i)
		sci_param->at(i) = _sci_param[i];

	s2p = RegisterObject<FRSIdParam>("FRS S2 ID Parameter");
	s3p = RegisterObject<FRSIdParam>("FRS S3 ID Parameter");
	sTar = RegisterObject<FRSTargetParam>("FRS S2 Target Parameter");

	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());
}

ClassImp(FRSIdParam);
ClassImp(FRSTargetParam);
ClassImp(RNFRSHit);
ClassImp(RNFRSHit::Id);

