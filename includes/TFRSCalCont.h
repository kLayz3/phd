#pragma once

#include "core/TContainer.hxx"
#include "nlohmann/json.hpp"
class TH2I;

struct RNSciCal {
	double E;
	std::vector<double> x;

	inline void Clean() noexcept { E = 0; x.clear(); }
	virtual ~RNSciCal() = default;
	ClassDef(RNSciCal, 1);
};

struct RNTPCCal {
	/* A single 'TPC measurement' means that
	 * some of the 4 anodes did a complete (x,y) measurement.
	 * Delay lines measure x - but we can associate them to anodes, 
	 * with a no-op. */
	
	struct Measurement {
		double x = NAN;
		double y = NAN;
		Measurement() = default;
		Measurement(double x, double y) : x(x), y(y) {}
		virtual ~Measurement() = default;
	};
	
	using Measurements = std::vector< 
		std::array<Measurement, 4>
	>;
	
	Measurements hits;
	inline void Clean() noexcept { hits.clear(); }
	virtual ~RNTPCCal() = default;
	ClassDef(RNTPCCal, 1);
};

struct alignas(util::CL) RNFRSCal {
	constexpr static i32 N_VALID_TPC = 4;

	std::array<RNTPCCal, N_VALID_TPC> tpc;
	RNSciCal sci21, sci22, sci31;

	inline void Clean() noexcept { 
		for(auto& t : tpc) t.Clean();
		sci21.Clean();
		sci22.Clean();
		sci31.Clean();
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

	template<std::size_t I>
    decltype(auto) get() &       noexcept { return get_helper<I>(*this); }
    template<std::size_t I> 
    decltype(auto) get() const & noexcept { return get_helper<I>(*this); }
    template<std::size_t I> 
    decltype(auto) get() &&      noexcept { return get_helper<I>(std::move(*this)); }
    template<std::size_t I>
    decltype(auto) get() const&& noexcept { return get_helper<I>(std::move(*this)); }

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

namespace std {
	template<> struct tuple_size<TPCParam> : integral_constant<size_t, TPCParam::N_PARAMS> {};
	template<> struct tuple_element<0, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<1, TPCParam> { using type = array<double, 2>; };
	template<> struct tuple_element<2, TPCParam> { using type = array<double, 4>; };
	template<> struct tuple_element<3, TPCParam> { using type = array<double, 4>; };
	template<> struct tuple_element<4, TPCParam> { using type = array<array<double, 2>, 4>; };
	template<> struct tuple_element<5, TPCParam> { using type = array<array<double, 2>, 4>; };
}

class TFRSCalCont: public TContainer<RNFRSCal> {
public:
	inline static nlohmann::json setup {}; 
	inline static std::array<TPCParam, 7> _tpc_param {};
		
	TH2I* h1_xy_s2_before_target;
	TH2I* h1_ab_s2_before_target;
	TH2I* h1_xy_s2_after_target;
	TH2I* h1_ab_s2_after_target;

	std::string *setupName{};
	std::array<TPCParam, 7> *tpc_param{};

	TFRSCalCont();
	void Setup() override;
	void Init(TDictInfo info) override;
};
