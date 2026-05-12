#pragma once

#include "TFOOTMapCont.h"
#include "monad/monad.hxx"
#include "util/json_struct_def.hh"
#include "TFOOTCalCont.h"

#include "util/HitMatrix.hxx"
#include "util/FTrack.h"

class TH2I;

struct FOOTHit {
	double Q; 
	double m; // measurement [mm]
	
	FOOTHit() = default;
	FOOTHit(double Q_, double m_) :
		Q(Q_), m(m_) {}

	virtual ~FOOTHit() = default;
	ClassDef(FOOTHit, 1);
};

struct RNFOOTPair {
	std::vector<FOOTHit> x;
	std::vector<FOOTHit> y;
	double z;
	RNFOOTPair() = default;

	inline void Clean() noexcept { x.clear(); y.clear(); }
	virtual ~RNFOOTPair() = default;
	ClassDef(RNFOOTPair, 1);
};
extern template struct HitMatrix<RNFOOTPair>; // instantiated in TFOOTHitProc.cxx
extern template struct Track<RNFOOTPair>; // instantiated in TFOOTHitProc.cxx

struct RNFOOTTrack {
	double x0, y0;
	double theta, phi;
	double Q;

	RNFOOTTrack() = default;
	virtual ~RNFOOTTrack() = default;
	ClassDef(RNFOOTTrack, 1);
};

struct RNFOOTHit {
	static constexpr u32 N_PAIRS = N_FOOT_DETECTORS / 2;
	std::array<RNFOOTPair, N_PAIRS> pair;
	std::vector<RNFOOTTrack> track;

	RNFOOTHit() = default;

	inline void Clean() noexcept { 
		for(auto& p: pair) p.Clean(); 
		track.clear(); 
	}
	virtual ~RNFOOTHit() = default;
	ClassDef(RNFOOTHit, 1);
};

struct TFOOTHitCont : TContainer<RNFOOTHit> {
	static constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	inline static nlohmann::json setup {};
	inline static FOOTBoxParam _box;
	
	FOOTBoxParam* box;
	std::string* setupName;
	TH2I* h_corr_all[N_PAIRS];
	TH2I* h_corr_all_sorted[N_PAIRS];
	TH2I* h_corr_gud[N_PAIRS];
	TH1I* h_single_all[N_FOOT_DETECTORS];
	TH1I* h_single_gud[N_FOOT_DETECTORS];
	TFOOTHitCont();

	void Setup() override;
	void Init(TDictInfo ) override;
};
