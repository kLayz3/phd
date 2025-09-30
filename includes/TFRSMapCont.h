#pragma once

#include "TContainer.hxx"

struct RNSciMap {
	static constexpr int MAX_SIZE = 10;
	std::array<std::vector<i32>, 2> tdc{}; // [0] = left, [1] = right;
	u16 qdc[2]{};                          // [0] = left, [1] = right;
	
	RNSciMap() {
		for(auto& t : tdc) t.reserve(MAX_SIZE);
	}
	inline void Clean() noexcept { 
		for(auto& t : tdc) t.clear(); 
		memset(qdc, 0, sizeof(qdc)); 
	}
	virtual ~RNSciMap() = default;
	ClassDef(RNSciMap, 1);
};

struct RNTPCMap {
	static constexpr int MAX_SIZE = 64;
	std::array<std::vector<i32>, 2> tdc_l{};
	std::array<std::vector<i32>, 2> tdc_r{};
	std::array<std::vector<i32>, 4> tdc_a{};
	std::vector<i32> tdc_ref{};
	u16 adc[4]{};

	RNTPCMap() { 
		for(auto& t : tdc_l) 
			t.reserve(MAX_SIZE);
		for(auto& t : tdc_r)
			t.reserve(MAX_SIZE);
		for(auto& t : tdc_a)
			t.reserve(MAX_SIZE);
		tdc_ref.reserve(MAX_SIZE);
	}
	inline void Clean() noexcept { 
		for(auto& t : tdc_l)
			t.clear();
		for(auto& t : tdc_r)
			t.clear();
		for(auto& t : tdc_a) 
			t.clear();
		tdc_ref.clear();

		memset(adc, 0, sizeof(adc));
	}
	virtual ~RNTPCMap() = default;
	ClassDef(RNTPCMap, 1);
};

template<uint32_t N>
struct RNMUSICMap {
	static_assert(N > 1 and N < 100, "Number of RNMUSICMap anodes >100 or <2 ?");
	static constexpr int size() { return static_cast<int>(N); }

	/* No multihit yet, TDC measurements also garbage. Just take ADC. */
	u16 e[N];
	RNMUSICMap() = default;
	inline void Clean() noexcept { memset(e, 0, sizeof(e)); }

	virtual ~RNMUSICMap() = default;
	ClassDef(RNMUSICMap, 1);
};

struct RNFRSMap {
	std::array<RNSciMap, 4> sci;
	std::array<RNTPCMap, 7> tpc;
	std::array<RNMUSICMap<8>, 2> music;
	uint32_t tpat;

	inline void Clean() noexcept {
		for(auto& s : sci) s.Clean();
		for(auto& t : tpc) t.Clean();
		for(auto& m : music) m.Clean();
		tpat = static_cast<u32>(-1);
	}
	virtual ~RNFRSMap() = default;
	ClassDef(RNFRSMap, 1);
};

class TFRSMapCont : public TContainer<RNFRSMap> {
	static_assert(sizeof(Int_t) == sizeof(i32), 
		"`i32` and `Int_t` unequal size? Change the std::memcpy to something human in the ProcessEntry!\n");
public: 
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
	
	TH1I* h1_sci_ml[4];
	TH1I* h1_sci_mr[4];
	TH1I* h1_tpc_ml[7];
	TH1I* h1_tpc_mr[7];
	TH1I* h1_tpc_ma1[7];
	TH1I* h1_tpc_ma2[7];
	TH1I* h1_tpc_csum[7][4];

	TFRSMapCont();
};
