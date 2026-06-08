#include "TFOOTHitCont.h"
#include "TFOOTHitProc.h"
#include "TH2I.h"
#include "TH1I.h"
#include "util/JSONParser.h"

using json = nlohmann::json;

std::ostream& operator<<(std::ostream& os, const FOOTHit& rhs) noexcept {
	os << mnd::msg("(%s%.2f%s,%d,%s%5.1f%s)", 
		KBH_CYN, rhs.Q.q, KNRM,
		rhs.Q.fCM,
		KBH_RED, rhs.m, KNRM);
	return os;
}

TFOOTHitCont::TFOOTHitCont() : TContainer("FOOT") {}

RNFOOTTrack::RNFOOTTrack (
	const std::array<double, 2>& xline, 
	const std::array<double, 2>& yline, 
	double Q_, 
	double score_,
	std::size_t n_
#ifdef MND_FOOTTRACK_DEBUG
		,
		const std::array<double, N_PAIRS>& _ox, 
		const std::array<double, N_PAIRS>& _oy, 
		const std::array<double, N_PAIRS>& _oz,
		const std::array<double, N_PAIRS>& _oq,
		const std::array<double, N_PAIRS>& _osq
#endif
) : x0(xline[0]), y0(yline[0]), ax(xline[1]), ay(yline[1]),
	Q(Q_), score(score_), n(n_) 
#ifdef MND_FOOTTRACK_DEBUG
		,
		_x(_ox), 
		_y(_oy), 
		_z(_oz),
		_q(_oq),
		_sq(_osq)
#endif

{}

void TFOOTHitCont::Init(TDictInfo info) {
auto it = info.find("Setup");	
	if(it == info.end()) 
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	
	setup = ParseJSON(file_name);
	setup["file_name"] = file_name;
	
	if(!setup.contains("box")) ERROR("\'box\' key not found in JSON file: %s\n", file_name.c_str());
		
	UNROLL_JSON_PARAM(_box, setup["box"], 7);
}

using FOOTParams = std::array<FOOTParam, 2>;
template<> void Add(FOOTParams&, const FOOTParams&) {}

void TFOOTHitCont::Setup() {
	h1_qtrack = RegisterObject<TH1I>("h1_qtrack", "Charge (Q) of recognized tracks",
		100, 0, 8);
	h1_track_nsampled = RegisterObject<TH1I>("h1_qtrack_nsampled", "Points sampled (N) in the track",
		12, -0.5, 5.5);
	
	setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());
	box = RegisterObject<FOOTBoxParam>("box", _box);

	/* TFOOTHitCont cannot know ahead of time which FOOT will get put into which position,
	 * for the tracking. It's the Proc which decides that.
	 * Therefore, we allocate the `foot_param` objects here, but leave them defaulted.
	 * The TFOOTHitProc ctor will fill the rest. */

	for(u32 i=0; i<N_PAIRS; ++i) {
		foot_param[i] = RegisterObject<FOOTParams>(Form("%u_setup", i), {});
	}
}

ClassImp(FOOTQ);
ClassImp(FOOTHit);
ClassImp(RNFOOTPair);
ClassImp(RNFOOTTrack);
ClassImp(RNFOOTHit);
