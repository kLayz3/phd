#include "TFRSCalCont.h"
#include <filesystem>
#include <regex>

#include "util/JSONParser.h"
#include "util/PolyFitter.h"

using nlohmann::json;
namespace fs = std::filesystem;

static nlohmann::json setup {};

bool SCIDEIntoQConverter::matches_file(std::string_view fname) const {
    /* `fname` could be with an extension, or with fullpath appended.
     * In this case, just strip it out. */
    const std::string stem = fs::path(fname).stem().string();
    std::regex re;
    try {
        re = std::regex{this->regex};
    } catch(const std::exception& e) {
        ERROR("SCIDEIntoQConverter::matches_file(...): "
            "Compiling underlying regex: \'%s\' failed. Info: %s\n",
            this->regex.c_str(), e.what());
    }

    /* Regex needs to match entirely on the stem. */
    return std::regex_match(fname.begin(), fname.end(), re);
}
double SCIDEIntoQConverter::Q(const RNSciCal& s) const noexcept {
    if(!this->is_initialized_)
        QParamInit();
    const f64 de_l = std::max( (s.El - pedestal.left),  BELOW_PEDESTAL_VAL);
    const f64 de_r = std::max( (s.Er - pedestal.right), BELOW_PEDESTAL_VAL);
    const f64 inv = std::sqrt( de_l * de_r );
    return std::pow(f_ * inv, c_);
}

void SCIDEIntoQConverter::QParamInit() const {
    if( !mnd::isfinite(pedestal.left, pedestal.right) ) {
        ERROR("SCI: de-to-q converter, pedestal left parameter is null.\n");
        std::cerr << pedestal << std::endl;
    }
        
    std::vector<f64> x, y;
    for(auto [Q, qdc_mean] : this->values) {
        if(qdc_mean <= pedestal.left || qdc_mean <= pedestal.right)
            ERROR("SCI: de-to-q converter, for charge %d, its mean value <= pedestal?\n", Q);
        
        x.push_back( std::log(Q) );
        y.push_back( std::log(qdc_mean) );
    }
    auto r = PolyFit<1>(x,y);
    this->f_ = std::exp(-r[0]);
    this->c_ = 1.0 / r[1];
    this->is_initialized_ = true;
}

double SCIParam::Q(const RNSciCal& s) const noexcept {
    if(!current_converter)
        return NAN;
    return current_converter->Q(s);
}

/* Regex and fs are called in a small loop, but it doesn't really matter.
 * This call should only be at the init, not in some kind of a loop. */
u32 SCIParam::SetConverter(std::string_view fname) const {
    u32 cnt{0};
    for(const auto& cvt : this->de_to_q) {
        if(cvt.matches_file(fname)) {
            current_converter = &cvt;
            ++cnt;
        }
    }
    return cnt;
}
SCIDEIntoQConverter const* SCIParam::GetConverter() const {
    return current_converter;
}

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
		if(RNFRSCal::tpc_moniker.find(_tpc_i) == RNFRSCal::tpc_moniker.end())
			ERROR("TPC parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_tpc_i.c_str(), file_name.c_str(), _tpc_param.size());
		
		u32 i = RNFRSCal::tpc_moniker.at(_tpc_i);
		if(i >= RNFRSCal::N_VALID_TPC) continue;
		UNROLL_JSON_PARAM(_tpc_param[i], params, 9)
	}

	for(const auto& [_sci_i, params] : setup.at("SCI").items()) {
		if(RNFRSCal::sci_moniker.find(_sci_i) == RNFRSCal::sci_moniker.end())
			ERROR("Sci parameter named \'%s\' found in the %s JSON parameter file isn't mapped to 0..%zu index.",
				_sci_i.c_str(), file_name.c_str(), _sci_param.size());
		
		u32 i = RNFRSCal::sci_moniker.at(_sci_i);
		if(i >= RNFRSCal::N_VALID_SCI) continue;
		UNROLL_JSON_PARAM(_sci_param[i], params, 4)
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
ClassImp(SCIQDCPedestal);
ClassImp(SCIMeanQDC);
ClassImp(SCIDEIntoQConverter);
ClassImp(SCIParam);
