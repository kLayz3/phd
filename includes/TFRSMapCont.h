#pragma once

#include "monad/monad.hxx"
#include "TH1I.h"
#include "TH2I.h"

struct RNSciMap {
	static constexpr int MAX_SIZE = 10;
	
	struct Measurement {
		static constexpr int TDC_INVALID = -1;
		int tdc_l = TDC_INVALID;
		int tdc_r = TDC_INVALID;
		
		Measurement() = default;
		virtual ~Measurement() = default;

		ClassDef(Measurement, 1);
	};

	std::vector<Measurement> tdc{}; 
	u16 qdc[2]{}; // [0] = left, [1] = right;
	
	RNSciMap() { tdc.reserve(MAX_SIZE); }
	inline void Clean() noexcept { 
		tdc.clear(); 
		memset(qdc, 0, sizeof(qdc)); 
	}

	virtual ~RNSciMap() = default;
	ClassDef(RNSciMap, 1);
};

/* The structs are layed out the following way to help cache locality.
 * Initially the layout was:
 *	std::array<std::vector<i32>, 2> tdc_l{};
 *	std::array<std::vector<i32>, 2> tdc_r{};
 *	std::array<std::vector<i32>, 4> tdc_a{};
 *	std::vector<i32> tdc_ref{};
 * -> This is diabolically bad, since we always take differences `tdc_l[??] - tdc_r[??]` and these references
 *  aren't on the same cache line. Worst yet - this heap allocated sequence could share close line with some FOOT
 *  data too, making congestion insane, and basically most values aren't (initially) cached in L1d.
 *  Having this layout and doing the simple loop above, slows down the whole program by a factor 100 !
 *
 *  A great read to the whole L1d locality + cache congestion and how indeed malloc works under the hood:
 *  https://people.freebsd.org/~lstewart/articles/cpumemory.pdf
 */

struct RNTPCMap {
	static constexpr i32 MAX_SIZE = 6;
	static_assert(MAX_SIZE > 0 && MAX_SIZE < 64, "64 is what Go4 gives us. Don't make higher capacity!");

	struct Measurement {
		static constexpr i32 TDC_INVALID = -1;
		i32 tdc_l = TDC_INVALID;     /* 4 bytes. */
		i32 tdc_r = TDC_INVALID;     /* 4 bytes. */
		std::array<i32, 2> tdc_a { TDC_INVALID, TDC_INVALID }; /* 8 bytes */ 

		Measurement() = default;
		virtual ~Measurement() = default;
		ClassDef(Measurement, 1);
	}; static_assert(sizeof(Measurement) == 24, "Huh?");

	std::array <
		std::vector<Measurement>, 2
	> tdc{};
	std::vector<i32> tdc_ref {};
	u16 adc[4] {};

	RNTPCMap() { for(auto& dl : tdc) dl.reserve(MAX_SIZE); tdc_ref.reserve(MAX_SIZE); }
	inline void Clean() noexcept {
		for(auto& tdc_dl : tdc) tdc_dl.clear();
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

struct alignas(mnd::CL) RNFRSMap {
	constexpr static i32 N_VALID_SCI = 4;
	constexpr static i32 N_VALID_TPC = 7;

	std::array<RNSciMap, N_VALID_SCI> sci;
	std::array<RNTPCMap, N_VALID_TPC> tpc;
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
	constexpr static i32 N_VALID_SCI = RNFRSMap::N_VALID_SCI;
	constexpr static i32 N_VALID_TPC = RNFRSMap::N_VALID_TPC;

	TH1I* h1_sci_ml[N_VALID_SCI];
	TH1I* h1_sci_mr[N_VALID_SCI];
	TH1I* h1_sci_diff_lr[N_VALID_SCI];

	TH1I* h1_tpc_araw[N_VALID_TPC][4];
	TH1I* h1_tpc_dl_lraw[N_VALID_TPC][2];
	TH1I* h1_tpc_dl_rraw[N_VALID_TPC][2];

	TH1I* h1_tpc_ml[N_VALID_TPC];
	TH1I* h1_tpc_mr[N_VALID_TPC];
	TH1I* h1_tpc_ma1[N_VALID_TPC];
	TH1I* h1_tpc_ma2[N_VALID_TPC];
	TH1I* h1_tpc_csum[N_VALID_TPC][4];
	TH1I* h1_tpc_ydiff[N_VALID_TPC][4];
	TH1I* h1_tpc_adiff[N_VALID_TPC][2]; /* Same TPC, anode (1)-(0) on same delay-line. */
	TH1I* h1_tpc_ldiff[N_VALID_TPC];    /* Same TPC, delay-left (1)-(0). */ 
	TH1I* h1_tpc_rdiff[N_VALID_TPC];    /* Same TPC, delay-left (1)-(0). */

	void Setup() override;

	TFRSMapCont();
};
