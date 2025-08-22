#pragma once
#include "TContainer.h"
#include "libs.hh"
#include "TH2I.h"
#include "TH2D.h"
#include <array>
#include "TOnce.hxx"

/* These defines are for loop-unrolling. */
#define _FOOT_N_STRIPS         640
#define _FOOT_N_ASIC            10
#define _FOOT_N_STRIPS_PER_ASIC 64

class TFOOTPedestalCont : public TContainer {
public:
	static constexpr int N_STRIPS = _FOOT_N_STRIPS; //!
	static constexpr int N_ASIC = _FOOT_N_ASIC; //!
	static constexpr int N_STRIPS_PER_ASIC = _FOOT_N_STRIPS_PER_ASIC; //!

public: // Inputs from Go4.
	u32* _FOOT; //!
	u32* _FOOTE; //! is array: [N_STRIPS]

	i32 FOOT_N; //!

public:
	// Objects held by the TContainer::_vc  
	TH2I* h2_raw;  //!
	TH2D* h2_mid;  //!
	TH2D* h2_corr; //!
	TGraph* gr_s0; //!
	TGraph* gr_s1; //!
	std::array<double, N_STRIPS> *gped, *gped_s, *gped_sf; //! indices [0], [1], [2], ..., [639]
	std::vector<int> *bad_strips; //! [0, 1, 2, ..., 639]
	TH2D* h2_ped_off_med; //!
	TH2D* h2_ped_off_avg; //!
	TH2D* h2_ped_off_diff; //!

public: // Exported back to ROOT

	/* Final ADC value after all three pedestal removals. */
	Double_t FOOTE[N_STRIPS];

/*
	// Category label for each strip in each event.
	enum SignalType {
		kBAD   = 0,
		kNOISE = 1,
		kHIT   = 2
	} stype[N_STRIPS];

	i32 n_bad = 0;
	i32 n_fired = 0;
	std::vector<int> fired_strips;
*/
public:

	TFOOTPedestalCont();
	TFOOTPedestalCont(int);
	virtual ~TFOOTPedestalCont();
	
	void SetId(int);
	
	void Clean(Option_t* option="") noexcept /* override */;
	void Init(TDictInfo info) /* override */;
	
	DECL_CONTAINER_METHODS

	ClassDef(TFOOTPedestalCont, 1);
};
