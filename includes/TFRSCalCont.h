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
		double x = NAN;
		double y = NAN;
		int ref_tdc = -1;

		Measurement() = default;
		Measurement(double x, double y, int r) : 
			x(x), y(y), ref_tdc(r) {}
		virtual ~Measurement() = default;
	}; /* ^^^ Per anode measurement. */
	
	using Measurements = std::vector <
		std::array<Measurement, 4>
	>;
	
	Measurements hits;
	inline void Clean() noexcept { hits.clear(); }
	virtual ~RNTPCCal() = default;
	ClassDef(RNTPCCal, 1);
};

struct alignas(mnd::CL) RNFRSCal {
	constexpr static i32 N_VALID_SCI = RNFRSMap::N_VALID_SCI;
	constexpr static i32 N_VALID_TPC = RNFRSMap::N_VALID_TPC;

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

struct TPCParam {
	constexpr static std::size_t N_PARAMS = 6;

	std::array<double, 2> x_offset;
	std::array<double, 2> x_factor;
	std::array<double, 4> y_offset;
	std::array<double, 4> y_factor;
	std::array<std::array<double, 2>, 4> csum_lim;
	std::array<std::array<double, 2>, 4> sci_ref_lim; // for y-calculation.

	GET_HELP_AUX_IMPL

	virtual ~TPCParam() = default;
	ClassDef(TPCParam, 1);

private:
	template<std::size_t I, typename Self>
	static decltype(auto) get_helper(Self&& self) noexcept {
		/* Brackets around return value because of decltype deduction rules. Check link:
		 * https://stackoverflow.com/questions/27557369/why-does-decltypeauto-return-a-reference-here */
		if constexpr(I == 0)      return (std::forward<Self>(self).x_offset);
		else if constexpr(I == 1) return (std::forward<Self>(self).x_factor);
		else if constexpr(I == 2) return (std::forward<Self>(self).y_offset);
		else if constexpr(I == 3) return (std::forward<Self>(self).y_factor);
		else if constexpr(I == 4) return (std::forward<Self>(self).csum_lim);
		else if constexpr(I == 5) return (std::forward<Self>(self).sci_ref_lim);
		else static_assert(I < N_PARAMS, "Index out of bounds.");
	} 
};

struct SCIParam {
	constexpr static std::size_t N_PARAMS = 3;
	constexpr static double channel_to_ns = 0.025;

	double x_offset;
	double x_factor;
	std::array<double, 2> cdiff_lim;

	GET_HELP_AUX_IMPL

	virtual ~SCIParam() = default;
	ClassDef(SCIParam, 1);

private:
	template<std::size_t I, typename Self>
	static decltype(auto) get_helper(Self&& self) noexcept {
		if constexpr(I == 0)      return (std::forward<Self>(self).x_offset);
		else if constexpr(I == 1) return (std::forward<Self>(self).x_factor);
		else if constexpr(I == 2) return (std::forward<Self>(self).cdiff_lim);
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
	template<> struct tuple_element<5, TPCParam> { using type = array<array<double, 2>, 4>; };

	template<> struct tuple_size<SCIParam> : integral_constant<size_t, SCIParam::N_PARAMS> {};
	template<> struct tuple_element<0, SCIParam> { using type = double; };
	template<> struct tuple_element<1, SCIParam> { using type = double; };
	template<> struct tuple_element<2, SCIParam> { using type = array<double, 2>; };
}

struct TFRSCalCont : TContainer<RNFRSCal> {

	inline static nlohmann::json setup {}; 
	inline static std::array<TPCParam, RNFRSCal::N_VALID_TPC> _tpc_param {};
	inline static std::array<SCIParam, RNFRSCal::N_VALID_SCI> _sci_param {};
	
	// Which name corresponds to which index in later naming convention.
	inline static const std::map<std::string, u32> tpc_moniker { 
		{"21", 0}, {"22", 1}, {"23", 2}, {"24", 3}, {"31", 4}, {"41", 5}, {"42", 6}
	}; 
	inline static const std::map<std::string, u32> sci_moniker { 
		{"21", 0}, {"22", 1},                       {"31", 2}, {"41", 3}
	};

	TH1I* h1_tpc_s2_before_target_nhit;
	TH2I* h2_xy_s2_before_target;
	TH2I* h2_ab_s2_before_target;

	TH1I* h1_tpc_s2_after_target_nhit;
	TH2I* h2_xy_s2_after_target;
	TH2I* h2_ab_s2_after_target;
	
	TH1I* h1_x_sc21_before_target;
	TH1I* h1_x_sc22_after_target;

	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_param{};
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_param{};
	std::string *setupName{};

	TFRSCalCont();
	void Setup() override;
	void Init(TDictInfo info) override;
};
