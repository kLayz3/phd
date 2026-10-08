#pragma once

#include "monad/monad.hxx"
#include "util/json_struct_def.hh"
#include "TFRSMapCont.h"

#include "TGraph.h"

class TH2I;

struct RNSciCal {
	struct Measurement {
		double x = NAN;
		double t = NAN;
		Measurement() = default;
		Measurement(double x, double t) : x(x), t(t) {}
		
		virtual ~Measurement() = default;
		ClassDef(Measurement, 1);
	};

    double El = NAN; // left  PMT QDC measurement
    double Er = NAN; // right PMT QDC measurement
	double E  = NAN; // quick n' dirty: sqrt(El * Er)
	std::vector<Measurement> hits;

	/* Return true if more than 1 hit registered by the MHTDC. */
	inline bool IsOk() const noexcept  { return hits.size() > 0; }
	inline void Clean() noexcept { E = (El = (Er = NAN)); hits.clear(); }
	virtual ~RNSciCal() = default;
	ClassDef(RNSciCal, 1);
};

struct RNTPCCal {
	/* A single 'TPC measurement' means that some of the 
	 * four anodes did a complete (x,y) measurement.
	 * Delay lines measure x - but we can associate them to anodes, 
	 * with a no-op. */
	
	struct Measurement {
		static constexpr i32 TDC_INVALID = RNTPCMap::Measurement::TDC_INVALID;
		double x  = NAN;
		std::array<double, 2> y = {NAN, NAN};
		i32 ref   = TDC_INVALID;
		u8 trusted = 0; /* Basically, if in map step all the list sizes match, and we trust this measurement. */

		Measurement() = default;
		Measurement(double _x, double _y0, double _y1 , i32 _ref, u8 _trusted) : 
			x(_x), y{_y0, _y1}, ref(_ref), trusted(_trusted) {}
		virtual ~Measurement() = default;

		inline double X() const noexcept { return x; }
		inline double Y() const noexcept {
			int N=0; double y0 = 0.0;
			if(!std::isnan(y[0])) { ++N; y0 += y[0]; } 
			if(!std::isnan(y[1])) { ++N; y0 += y[1]; }
			return y0 / N;
		}
		inline int AnodeMask() const noexcept {
			int r = 0;
			if(!std::isnan(y[0])) r |= 0x1;
			if(!std::isnan(y[1])) r |= 0x2;
			return r;
		}
	}; /* ^^^ Per delay-line measurement. */

	/* TPC is represented as two independent xy-measurements coming from 2 delay lines. */
	using Measurements = std::array <
		std::vector<Measurement>, 2
	>;
	Measurements hits;
	
	/* Returns X-position of the first registered hit. Quiet NaN if there are no hits. */
	inline double X0() const noexcept {
		int N=0; double x0 = 0.0;
		if(hits[0].size() > 0) x0 += hits[0][0].X(), ++N;
		if(hits[1].size() > 0) x0 += hits[1][0].X(), ++N;
		return x0 / N;
	}
	/* Returns Y-position of the first registered hit. Quiet NaN if there are no hits. */
	inline double Y0() const noexcept {
		i32 w0 = 0, w1 = 0; 
		double y0 = 0, y1 = 0;
		if(hits[0].size() > 0) w0 = 1 + (hits[0][0].AnodeMask() == 0x3), y0 = hits[0][0].Y();
		if(hits[1].size() > 0) w1 = 1 + (hits[1][0].AnodeMask() == 0x3), y1 = hits[1][0].Y();
		return (w0*y0 + w1*y1) / (w0 + w1); 
	}

	inline void Clean() noexcept { for(auto& d : hits) d.clear(); }
	virtual ~RNTPCCal() = default;
	ClassDef(RNTPCCal, 1);
};

struct RNFRSCal {
	constexpr static i32 N_VALID_SCI = RNFRSMap::N_VALID_SCI;
	constexpr static i32 N_VALID_TPC = RNFRSMap::N_VALID_TPC;
    constexpr static double S2_LENGTH = 4560.0;

	// Which name corresponds to which index in later naming convention.
	// Note, we keep this to match Map step.
    static constexpr u32 SCI21_I = 0;
    static constexpr u32 SCI22_I = 1;
    static constexpr u32 SCI31_I = 2;
    static constexpr u32 SCI41_I = 3;

    static constexpr u32 TPC21_I = 0;
    static constexpr u32 TPC22_I = 1;
    static constexpr u32 TPC23_I = 2;
    static constexpr u32 TPC24_I = 3;
    static constexpr u32 TPC41_I = 4;
    static constexpr u32 TPC42_I = 5;
    static constexpr u32 TPC31_I = 6;

	inline static const std::map<std::string, u32> tpc_moniker {
		{"21", TPC21_I},
        {"22", TPC22_I},
        {"23", TPC23_I},
        {"24", TPC24_I},
        {"41", TPC41_I},
        {"42", TPC42_I},
        {"31", TPC31_I}
	};
	static constexpr std::array<const char*, N_VALID_TPC> tpc_label = {
		"21", "22", "23", "24", "41", "42", "31"
	};
    inline static const std::map<std::string, u32> sci_moniker {
        {"21", SCI21_I},
        {"22", SCI22_I},
        {"31", SCI31_I},
        {"41", SCI41_I}
    };
	static constexpr std::array<const char*, N_VALID_SCI> sci_label = {
        "21", "22", "31", "41"
    };

	std::array<RNSciCal, N_VALID_SCI> sci;
	std::array<RNTPCCal, N_VALID_TPC> tpc;
    RNTrigMap trig; // Just gets directly mapped from the map step.

	inline void Clean() noexcept { 
		for(auto& s : sci) s.Clean();
		for(auto& t : tpc) t.Clean();
	}
	virtual ~RNFRSCal() = default;
	ClassDef(RNFRSCal, 2);
};

struct TPCParam {
	using arr2 = std::array<double,2>;
	using arr4 = std::array<double,4>;
	using arr22 = std::array<std::array<double,2>, 2>;
	using arr24 = std::array<std::array<double,2>, 4>;

	GET_HELP_AUX_IMPL;

    /* Nominal value... not all that important. It's just for a more 
     * pedantic (x,y,z) measurement of a hit of an anode. */
    static constexpr double TPC_WIDTH = 70.0;
    static constexpr u32 N_UPSTREAM_TPC   = 3;
    static constexpr u32 N_DOWNSTREAM_TPC = 1;
    static constexpr u32 N_S2_TPC = N_UPSTREAM_TPC + N_DOWNSTREAM_TPC;

	ADD_SERIALIZABLE_FIELD(arr2,   x_offset,             {}, 0);
	ADD_SERIALIZABLE_FIELD(arr2,   x_factor,             {}, 1);
	ADD_SERIALIZABLE_FIELD(arr4,   y_offset,             {}, 2);
	ADD_SERIALIZABLE_FIELD(arr4,   y_factor,             {}, 3);
	ADD_SERIALIZABLE_FIELD(arr24,  csum_lim,             {}, 4);
	ADD_SERIALIZABLE_FIELD(arr22,  anode_diff_lim,       {}, 5);
	ADD_SERIALIZABLE_FIELD(arr2,   dl_left_diff_lim,     {}, 6);
	ADD_SERIALIZABLE_FIELD(arr2,   dl_right_diff_lim,    {}, 7);
	ADD_SERIALIZABLE_FIELD(arr24,  sci_ref_lim,          {}, 8);
	ADD_SERIALIZABLE_FIELD(double, z0,                    0, 9); 
	
    arr2 zDL() const {
		return {z0 - TPC_WIDTH/4,
		        z0 + TPC_WIDTH/4 };
    }

	virtual ~TPCParam() = default;
	ClassDef(TPCParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(TPCParam, 9);

/* This type is responsible for β-correction of the
 * dE measurement of scintillators. Hence, only will be utilized in
 * the next, hit step. */
struct SCIPrimary {
    GET_HELP_AUX_IMPL
	using A2 = std::array<double, 2>;

	// Nullable according to the `beta` field. Passing in beta outside of <0,1> in parameter file
	// is undefined behaviour.
    ADD_SERIALIZABLE_FIELD(f64, beta, NAN, 0); // Referent velocity at maximal matter
	ADD_SERIALIZABLE_FIELD(A2,  fit,  {},  1); // Fit for E vs. 1/β^2

	/* SCI's give the 1st calibrated QDC value as:
	 * sqrt((L - pedL) * (R - pedR)) .
	 * Calculate the correction factor of this energy based on the incoming particle β,
	 * if the β can be calculated. Passing in β<0 or a NAN, returns back a factor 1.0. */
	[[nodiscard]] double BetaCorrection(double b) const noexcept;
	bool IsOk() const noexcept;
	void Reset() const noexcept;

	[[nodiscard]] TGraph* GetGraph(
		int ndiv  = 100,
		double lo = 0.56,
		double hi = 0.95
	) const;

    SCIPrimary() = default;
	virtual ~SCIPrimary() = default;
	ClassDef(SCIPrimary, 1);

private:
	void Init() const noexcept;
	mutable mnd::cache<double> qdc_ref{}; //!
	static std::mutex mtx_; //!
};

ADD_JSON_TYPE_RESOLUTION(SCIPrimary, 1)

struct SCIQDCPedestal {
    GET_HELP_AUX_IMPL
	constexpr static f64 INVALID_PEDESTAL = NAN;

    ADD_SERIALIZABLE_FIELD(f64, left,  NAN, 0);
    ADD_SERIALIZABLE_FIELD(f64, right, NAN, 1);
	
	constexpr inline bool IsDefaulted() const noexcept {
		return !mnd::isfinite(left, right)
			|| left  == INVALID_PEDESTAL
			|| right == INVALID_PEDESTAL;
	}
	constexpr inline bool IsValid() const noexcept {
		return !IsDefaulted();
	}

    SCIQDCPedestal() = default;
	virtual ~SCIQDCPedestal() = default;
	ClassDef(SCIQDCPedestal, 1);
};
ADD_JSON_TYPE_RESOLUTION(SCIQDCPedestal, 1)
bool operator==(SCIQDCPedestal const& , SCIQDCPedestal const& );
bool operator!=(SCIQDCPedestal const& , SCIQDCPedestal const& );

/* We don't globally have only a single parameter for dE -> Q conversion,
 * There can be different QDC gains set for different files. We tag the sequence of files
 * with a string regex. */
struct SCIMeanQDC {
	GET_HELP_AUX_IMPL
    ADD_SERIALIZABLE_FIELD(i32, Q,        -1,  0);
    ADD_SERIALIZABLE_FIELD(f64, qdc_mean, NAN, 1);

    SCIMeanQDC() = default;
	virtual ~SCIMeanQDC() = default;
	ClassDef(SCIMeanQDC, 1);
};
ADD_JSON_TYPE_RESOLUTION(SCIMeanQDC, 1)
bool operator==(SCIMeanQDC const& , SCIMeanQDC const& );
bool operator!=(SCIMeanQDC const& , SCIMeanQDC const& );

struct SCIDEIntoQConverter {
	GET_HELP_AUX_IMPL
    constexpr static f64 BELOW_PEDESTAL_VAL = 0.66; // Some random small number

    using SCIMeanQDCSeq = std::vector<SCIMeanQDC>;
    ADD_SERIALIZABLE_FIELD(SCIQDCPedestal, pedestal, {},  0);
    ADD_SERIALIZABLE_FIELD(SCIMeanQDCSeq,  values,   {},  1);
 
    /* Main method: convert SCI energy (E) to nominal charge (Q)
	 * "Energy" is defined as E0 = sqrt( (E(l) - Ped(l)) * (E(r) - Ped(r)) ),
	 * times the beta-correction factor (~0.95 - 1.05).
	 * If the dependence is: E(Q) = A * Q^a,
     * then: f = 1/A, c = 1/a <=> Q(E) = (f * E)^c */
	double Q(double ) const noexcept;

	/* Return just the initial energy:
	 * E0 = sqrt( (E(l) - Ped(l)) * (E(r) - Ped(r)) ) */
	double E0(const RNSciCal& ) const noexcept;

	/* Sci left-measured energy (ADC units): E(l) - Ped(l). */
	double E0_l(const RNSciCal& ) const noexcept;

	/* Sci right-measured energy (ADC units): E(l) - Ped(l). */
	double E0_r(const RNSciCal& ) const noexcept;

	inline void Reset() const noexcept { fit_.reset(); }
	void Init(bool verbose = false) const;

	constexpr inline bool IsDefaulted() const noexcept {
		return pedestal.IsDefaulted();
	}
	constexpr inline bool IsValid() const noexcept {
		return !IsDefaulted();
	}

	[[ nodiscard ]] std::pair<TGraph*, TGraph*> GetGraph(
		int ndiv  = 500,
		double lo = 10,
		double hi = 4000
	) const;

    SCIDEIntoQConverter() = default;

protected:
	/* Some cached values for quick Q- calculation.
	 * NB: if the object is re-evaluted, the values *need* to be recomputed, but the default
	 * JSON propagator cannot know this. Meaning that `Reset()` has to be called manually. */
	struct fit_coeff {
		double f, c;
		fit_coeff() : f{NAN}, c{NAN} {};
	};
	mutable mnd::cache<fit_coeff> fit_; //!
	static std::mutex mtx_; //!

public:
	virtual ~SCIDEIntoQConverter() = default;
	ClassDef(SCIDEIntoQConverter, 1);
};
ADD_JSON_TYPE_RESOLUTION(SCIDEIntoQConverter, 1)

bool operator==(SCIDEIntoQConverter const& , SCIDEIntoQConverter const& );
bool operator!=(SCIDEIntoQConverter const& , SCIDEIntoQConverter const& );

struct SCIParam {
	GET_HELP_AUX_IMPL
	constexpr static double channel_to_ns = 0.025;
	using arr2 = std::array<double,2>;

	ADD_SERIALIZABLE_FIELD(double,                x_offset,  0,  0);
	ADD_SERIALIZABLE_FIELD(double,                x_factor,  1,  1);
	ADD_SERIALIZABLE_FIELD(arr2,                  cdiff_lim, {}, 2);
	ADD_SERIALIZABLE_FIELD(double,                z0,        0,  3);
	ADD_SERIALIZABLE_FIELD(SCIPrimary,            primary,   {}, 4);
    ADD_SERIALIZABLE_FIELD(SCIDEIntoQConverter,   qdc,       {}, 5);

    SCIParam() = default;

	/* Returns the SCI final calibrated charge (after pedestal subtraction and β-correction). */
    double Q(const RNSciCal& s, double beta = NAN) const noexcept;
	
	/* Returns the SCI final calibrated charge (after pedestal subtraction and β-correction).
	 * Input is the preliminary energy measurement: sqrt((E_l - ped_l) * (E_r - ped_r)) .*/
    double Q(double E0, double beta = NAN) const noexcept;

	/* Sci measured energy (ADC units) after pedestal subtraction. */
	double E0(const RNSciCal& ) const noexcept;

	/* Sci left-measured energy (ADC units) after pedestal subtraction. */
	double E0_l(const RNSciCal& ) const noexcept;

	/* Sci right-measured energy (ADC units) after pedestal subtraction. */
	double E0_r(const RNSciCal& ) const noexcept;

	/* Sci measured energy (ADC units) after pedestal subtraction and β-correction. */
	double E(const RNSciCal& , double beta=NAN ) const noexcept;

public:
	virtual ~SCIParam() = default;
	ClassDef(SCIParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(SCIParam, 5)

struct TrigParamSingle {
    GET_HELP_AUX_IMPL
    ADD_SERIALIZABLE_FIELD(std::string, label,            "",    0);
    ADD_SERIALIZABLE_FIELD(u32,         tpat,             0,     1);
    ADD_SERIALIZABLE_FIELD(u32,         ttype,            1,     2);
    ADD_SERIALIZABLE_FIELD(bool,        enabled,          true,  3);
    ADD_SERIALIZABLE_FIELD(bool,        contains_tracker, false, 4);

    bool IsActive(u32 ) const noexcept;
    TrigParamSingle() = default;
	virtual ~TrigParamSingle() = default;
	ClassDef(TrigParamSingle, 1);
};
ADD_JSON_TYPE_RESOLUTION(TrigParamSingle, 4)

struct TrigParam {
    GET_HELP_AUX_IMPL

    static constexpr u32 MAX_TPAT = 16;
    ADD_SERIALIZABLE_FIELD(std::vector<TrigParamSingle>, trigger_map, {},  0);

    bool HasFOOT(u32 ) const noexcept;

    using TLabels = std::array<std::string, MAX_TPAT+1>; // Rather return owned string to not worry about dangles.
    using TLabelsMap = std::map<u32, std::string>;
    TLabels GetTrigLabels() const;
    TLabelsMap GetTrigLabelsMap() const;

    TrigParam() = default;
	virtual ~TrigParam() = default;
	ClassDef(TrigParam, 1);
};
ADD_JSON_TYPE_RESOLUTION(TrigParam, 0)

inline void Add(TPCParam&, const TPCParam&) {}
inline void Add(SCIParam&, const SCIParam&) {}
inline void Add(TrigParam&, const TrigParam&) {}

struct TFRSCalCont : TContainer<RNFRSCal> {
	friend struct TFRSCalProc;
	friend struct TFRSHitProc;
	
	TH2I* h2_tpc_xy[RNFRSCal::N_VALID_TPC][2];
	TH1I* h1_tpc_y[RNFRSCal::N_VALID_TPC][4];
	TH1I* h1_tpc_mask[RNFRSCal::N_VALID_TPC][2];
	TH1I* h1_x_sc21_before_target;
	TH1I* h1_x_sc22_after_target;
	TH1I* h1_s3_tdc_dt;

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param{};
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_param{};
    TrigParam *trig_param{};
	std::string *setupName{};

	TFRSCalCont();
	void Setup() override;
	void Init(TDictInfo info) override;

    static std::array<double, TPCParam::N_S2_TPC> z_s2_tpc (
        std::array<TPCParam, RNFRSCal::N_VALID_TPC> *
    );
    static std::array<std::array<double, 2>, TPCParam::N_S2_TPC> z_s2_tpc_delay_lines(
        std::array<TPCParam, RNFRSCal::N_VALID_TPC> *
    );

private:
	static nlohmann::json setup;
};
