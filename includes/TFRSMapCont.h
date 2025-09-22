#pragma once

#include "TContainer.h"

class TFRSMapCont : public TContainer {
	static_assert(sizeof(Int_t) == sizeof(i32), 
		"`i32` and `Int_t` unequal size? Change the std::memcpy to something human in the ProcessEntry!\n");
public: 
	class Sci {
	public:
		static constexpr int MAX_SIZE = 10;
		Int_t* _nhit_raw[2]; //! [0] = left, [1] = right;
		Int_t* _data_raw[2]; //! [0] = left, [1] = right;
		Int_t* _qdc_raw[2];  //! [0] = left, [1] = right;
	public:
		std::array<std::vector<i32>, 2> tdc{}; // [0] = left, [1] = right;
		u16 qdc[2]{};                          // [0] = left, [1] = right;
		inline void __clear() { 
			for(int i=0; i<2; ++i) tdc[i].clear(); 
			memset(qdc, 0, sizeof(qdc)); 
		}
		
		Sci() {
			for(int i=0; i<2; ++i)
				tdc[i].reserve(MAX_SIZE);
		}
		virtual ~Sci() {};
		ClassDef(Sci, 1);
	};

	// --------------------------------------------------- //
	class TPC {
	public:
		static constexpr int MAX_SIZE = 64;
		Int_t *_tpc_aa{}; //!
		Int_t *_tpc_lt[2], *_tpc_rt[2], *_tpc_at[4]; //!
		Int_t *_tpc_ltn[2], *_tpc_rtn[2], *_tpc_atn[4]; //!
	public:
		static u32 tpc_ref_i[7]; //!

		std::array<std::vector<i32>, 2> tdc_l{};
		std::array<std::vector<i32>, 2> tdc_r{};
		std::array<std::vector<i32>, 4> tdc_a{};
		u16 adc[4]{};

		inline void __clear() { 
			for(int i=0; i<2; ++i) {
				tdc_l[i].clear();
				tdc_r[i].clear();
			}
			for(int i=0; i<4; ++i)
				tdc_a[i].clear();

			memset(adc, 0, sizeof(adc));
		}
		TPC() { 
			for(int i=0; i<2; ++i) {
				tdc_l[i].reserve(MAX_SIZE);
				tdc_r[i].reserve(MAX_SIZE);
			}
			for(int i=0; i<4; ++i)
				tdc_a[i].reserve(MAX_SIZE);
		}
		virtual ~TPC() {};
		ClassDef(TPC, 1);
	};

	// --------------------------------------------------- //

	template<u32 N>
	class MUSIC {
		static_assert(N > 1 and N < 100, "Number of MUSIC anodes >100 or <2 ?");
	public:
		/* No multihit yet, TDC measurements also garbage. Just take ADC. */
		Int_t* _music_raw; //!
		u16 e[N];
		inline void __clear() { memset(e, 0, sizeof(e)); }
		MUSIC() = default;
		virtual ~MUSIC() {};

		ClassDef(MUSIC,1);
	};

	// --------------------------------------------------- //
	
	Int_t* _pattern; //!
	uint32_t tpat;

	std::array<Sci, 4> sci;
	std::array<TPC, 7> tpc{};
	std::array<MUSIC<8>, 2> music{};
	
	TH1I* h1_sci_ml[4];
	TH1I* h1_sci_mr[4];
	TH1I* h1_tpc_ml[7]; 
	TH1I* h1_tpc_mr[7];
	TH1I* h1_tpc_ma1[7];
	TH1I* h1_tpc_ma2[7];

	void Clear();
	
	// Add histos I cba for now ...
public:
	TFRSMapCont();
	virtual ~TFRSMapCont();

	DECL_CONTAINER_METHODS

	ClassDef(TFRSMapCont, 1);
};
