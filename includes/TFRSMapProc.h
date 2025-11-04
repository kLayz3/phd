#pragma once

#include "core/TProcessor.hxx"
#include "TFRSMapCont.h"
#include "TFRSGo4Cont.hxx"

struct TFRSMapProc : TProcessor <
	TFRSMapCont
	(TFRSGo4Cont)
> {
	using Base = TProcessor<TFRSMapCont(TFRSGo4Cont)>;
	enum class DoAnalysis : u8 { YES, NO } do_analysis;

	TFRSMapProc() = default;
	TFRSMapProc(TFRSMapCont& out, const TFRSGo4Cont& in, DoAnalysis y_or_n) 
		: Base(out, in), do_analysis(y_or_n) {}
	
	void  ProcessEntry() noexcept; 
	void _ProcessEntry() noexcept; 

private:
	struct Sci {
		static constexpr int MAX_SIZE = RNSciMap::MAX_SIZE;
		Int_t* _nhit_raw[2]; // [0] = left, [1] = right;
		Int_t* _data_raw[2]; // [0] = left, [1] = right;
		Int_t* _qdc_raw[2];  // [0] = left, [1] = right;
	};

	// --------------------------------------------------- //
	struct TPC {
		static constexpr int MAX_SIZE = RNTPCMap::MAX_SIZE;
		static constexpr u32 tpc_ref[MAX_SIZE] {
			0, // TPC21 => SCI21
			0, // TPC22 => SCI21
			0, // TPC23 => SCI21
			1, // TPC24 => SCI22
			2, // TPC31 => SCI31
			3, // TPC41 => SCI41
			3  // TPC42 => SCI41
		};
		Int_t *_tpc_aa{};
		Int_t *_tpc_lt[2], *_tpc_rt[2], *_tpc_at[4];
		Int_t *_tpc_ltn[2], *_tpc_rtn[2], *_tpc_atn[4];
		Int_t *_sci_timerefn, *_sci_timeref;
	};

	struct MUSIC {
		Int_t* _music_raw;
	};

	// --------------------------------------------------- //
	
	Int_t* _pattern;	
	std::array<Sci, 4> sci;
	std::array<TPC, 7> tpc;
	std::array<MUSIC, 2> music;

	void SetupPointers();

	/* Some TPC vars' local temporary storage. */
	Int_t _nhits_l[2] {};
	Int_t _nhits_r[2] {};
	Int_t _nhits_a[4] {};
	Int_t _nhits_s {};

	std::array<Int_t, RNTPCMap::MAX_SIZE> _temp {};
};
