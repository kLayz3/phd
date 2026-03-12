#include "TFRSHitCont.h"
#include "util/JSONParser.h"

using nlohmann::json;

TFRSHitCont::TFRSHitCont() : TContainer("FRS") {}  
 
void TFRSHitCont::Init(TDictInfo info) {
	auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;

	setup = ParseJSON(file_name);
	setup["file_name"] = file_name;

	auto j_it = setup.find("FRS");
	if(j_it == setup.end()) ERROR("\'FRS\' key not found in JSON file: %s\n", file_name.c_str());
	auto& jfrs = j_it.value();

	j_it  = jfrs.find("S2");
	if(j_it == jfrs.end()) ERROR("\'S2\' key not found in \'FRS\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_s2p, j_it.value(), 4);
	
	/* Unpack the S2 target param values. */
	auto& js2 = j_it.value();
	j_it = j_it.value().find("target");
	if(j_it == js2.end()) ERROR("\'targte\' key not found in \'S2\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_sTar, j_it.value(), 1);
	
	j_it  = jfrs.find("S3");
	if(j_it == jfrs.end()) ERROR("\'S3\' key not found in \'FRS\' section of JSON file: %s\n", file_name.c_str());
	UNROLL_JSON_PARAM(_s3p, j_it.value(), 4);
}

void Add(FRSIdParam&, const FRSIdParam&) {}
void Add(FRSTargetParam&, const FRSTargetParam&) {}

void TFRSHitCont::Setup() {
	h2_track_x = RegisterObject<TH2D>("h2_s2_track_x", "S2 Tracking;z[mm];x[mm]", 400, 0, 4200, 200, -100, 100);
	h2_track_y = RegisterObject<TH2D>("h2_s2_track_y", "S2 Tracking;z[mm];y[mm]", 400, 0, 4200, 200, -100, 100);
	s2p = RegisterObject<FRSIdParam>("FRS S2 ID Parameter");
	s3p = RegisterObject<FRSIdParam>("FRS S3 ID Parameter");
	sTar = RegisterObject<FRSTargetParam>("FRS S2 Target Parameter");
	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());
}

ClassImp(FRSIdParam);
ClassImp(FRSTargetParam);
ClassImp(RNFRSHit);
ClassImp(RNFRSHit::Id);

