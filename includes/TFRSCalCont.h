#pragma once

#include "monad/monad.hxx"
#include "nlohmann/json.hpp"
#include "TFRSMapCont.h"

class TH2I;

#define GET_HELP_AUX_IMPL  \
	template<std::size_t I> \
    decltype(auto) get() &       noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() const & noexcept { return get_helper<I>(*this); } \
    template<std::size_t I>  \
    decltype(auto) get() &&      noexcept { return get_helper<I>(std::move(*this)); } \
    template<std::size_t I> \
    decltype(auto) get() const&& noexcept { return get_helper<I>(std::move(*this)); }

struct RNSciCal {
	struct Measurement {
		double x = NAN;
		double t = NAN;
		Measurement() = default;
		Measurement(double x, double t) : x(x), t(t) {}
		
		virtual ~Measurement() = default;
		ClassDef(Measurement, 1);
	};

	double E = NAN;
	std::vector<Measurement> hits;

	inline bool IsOk() const noexcept  { return hits.size() > 0; }
	inline void Clean() noexcept { E = 0; hits.clear(); }
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

struct alignas(mnd::CL) RNFRSCal {
	constexpr static i32 N_VALID_SCI = RNFRSMap::N_VALID_SCI;
	constexpr static i32 N_VALID_TPC = RNFRSMap::N_VALID_TPC;
	
	inline static constexpr std::array<const char*, N_VALID_TPC> tpc_label = {
		"21", "22", "23", "24", "41", "42", "31"
	};	
	inline static constexpr std::array<const char*, N_VALID_SCI> sci_label = {
		"21", "22",             "31", "41"
	};

	std::array<RNSciCal, N_VALID_SCI> sci;
	std::array<RNTPCCal, N_VALID_TPC> tpc;

	inline void Clean() noexcept { 
		for(auto& s : sci) s.Clean();
		for(auto& t : tpc) t.Clean();
	}
	virtual ~RNFRSCal() = default;
	ClassDef(RNFRSCal, 1);
};

/* As a sanity check, we put setup parameters in the container, 
 * not in the processor,
 * to have it as an object in the ROOT file for later. */

template<typename U, size_t N, size_t M>
using arr2d = std::array<std::array<U, N>, M>;

struct TPCParam {
	constexpr static std::size_t N_PARAMS = 10;

	std::array<double, 2> x_offset;
	std::array<double, 2> x_factor;
	std::array<double, 4> y_offset;
	std::array<double, 4> y_factor;
	std::array<std::array<double, 2>, 4> csum_lim;
	std::array<std::array<double, 2>, 2> anode_diff_lim;
	std::array<double, 2> dl_left_diff_lim;
	std::array<double, 2> dl_right_diff_lim;
	std::array<std::array<double, 2>, 4> sci_ref_lim; // for y-calculation.
	
	double z0; // nominal z position.
	
	GET_HELP_AUX_IMPL

	virtual ~TPCParam() = default;
	ClassDef(TPCParam, 1);

private:
	template<std::size_t I, typename Self>
	static decltype(auto) get_helper(Self&& self) noexcept {
		/* Brackets around return value because of decltype deduction rules. Check link:
		 * https://stackoverflow.com/questions/27557369/why-does-decltypeauto-return-a-reference-here */
		if      constexpr(I == 0) return (std::forward<Self>(self).x_offset);
		else if constexpr(I == 1) return (std::forward<Self>(self).x_factor);
		else if constexpr(I == 2) return (std::forward<Self>(self).y_offset);
		else if constexpr(I == 3) return (std::forward<Self>(self).y_factor);
		else if constexpr(I == 4) return (std::forward<Self>(self).csum_lim);
		else if constexpr(I == 5) return (std::forward<Self>(self).anode_diff_lim);
		else if constexpr(I == 6) return (std::forward<Self>(self).dl_left_diff_lim);
		else if constexpr(I == 7) return (std::forward<Self>(self).dl_right_diff_lim);
		else if constexpr(I == 8) return (std::forward<Self>(self).sci_ref_lim);
		else if constexpr(I == 9) return (std::forward<Self>(self).z0);
		else static_assert(I < N_PARAMS, "Index out of bounds.");
	} 
};

struct SCIParam {
	constexpr static std::size_t N_PARAMS = 4;
	constexpr static double channel_to_ns = 0.025;

	double x_offset;
	double x_factor;
	std::array<double, 2> cdiff_lim;
	
	double z0; // nominal z position, w.r.t. 'Abstaende S2 Strahlzeit 2024'.
	GET_HELP_AUX_IMPL

	virtual ~SCIParam() = default;
	ClassDef(SCIParam, 1);

private:
	template<std::size_t I, typename Self>
	static decltype(auto) get_helper(Self&& self) noexcept {
		if      constexpr(I == 0) return (std::forward<Self>(self).x_offset);
		else if constexpr(I == 1) return (std::forward<Self>(self).x_factor);
		else if constexpr(I == 2) return (std::forward<Self>(self).cdiff_lim);
		else if constexpr(I == 3) return (std::forward<Self>(self).z0);
		else static_assert(I < N_PARAMS, "Index out of bounds.");
	}
};

namespace std {
	template<> struct tuple_size<TPCParam> : integral_constant<size_t, TPCParam::N_PARAMS> {};
	template<> struct tuple_element<0, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<1, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<2, TPCParam> { using type = array<double, 4>; };
	template<> struct tuple_element<3, TPCParam> { using type = array<double, 4>; };
	template<> struct tuple_element<4, TPCParam> { using type = array<array<double, 2>, 4>; };
	template<> struct tuple_element<5, TPCParam> { using type = array<array<double, 2>, 2>; };
	template<> struct tuple_element<6, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<7, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<8, TPCParam> { using type = array<array<double, 2>, 4>; };
	template<> struct tuple_element<9, TPCParam> { using type = double; };

	template<> struct tuple_size<SCIParam> : integral_constant<size_t, SCIParam::N_PARAMS> {};
	template<> struct tuple_element<0, SCIParam> { using type = double; };
	template<> struct tuple_element<1, SCIParam> { using type = double; };
	template<> struct tuple_element<2, SCIParam> { using type = array<double, 2>; };
	template<> struct tuple_element<3, SCIParam> { using type = double; };
}

struct TFRSCalCont : TContainer<RNFRSCal> {

	inline static nlohmann::json setup {}; 
	inline static std::array<TPCParam, RNFRSCal::N_VALID_TPC> _tpc_param {};
	inline static std::array<SCIParam, RNFRSCal::N_VALID_SCI> _sci_param {};
	
	// Which name corresponds to which index in later naming convention.
	// Note, we keep this to match Go4.
	inline static const std::map<std::string, u32> tpc_moniker { 
		{"21", 0}, {"22", 1}, {"23", 2}, {"24", 3}, {"41", 4}, {"42", 5}, {"31", 6}
	}; 
	inline static const std::map<std::string, u32> sci_moniker { 
		{"21", 0}, {"22", 1},                       {"31", 2}, {"41", 3}
	};

	TH2I* h2_tpc_xy[RNFRSCal::N_VALID_TPC][2];
	TH1I* h1_tpc_y[RNFRSCal::N_VALID_TPC][4];

	TH1I* h1_tpc_mask[RNFRSCal::N_VALID_TPC][2];

	TH1I* h1_x_sc21_before_target;
	TH1I* h1_x_sc22_after_target;

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param{};
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_param{};
	std::string *setupName{};

	TFRSCalCont();
	void Setup() override;
	void Init(TDictInfo info) override;
};
