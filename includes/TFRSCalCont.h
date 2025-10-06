#pragma once

#include "TContainer.hxx"
#include "nlohmann/json.hpp"
#include <array>
class TH2I;

struct alignas(util::CL) RNFRSCal {
	class Position {
	public:
		double x{}, y{};
		double a{}, b{};
		Position() = default;
		virtual ~Position() = default;
		ClassDef(Position, 1);
	};
	using Positions = std::vector<Position>;
	
	Positions before_target{};
	Position after_target{};
	Position s3 {};

	RNFRSCal() = default;
	virtual ~RNFRSCal() = default;
	ClassDef(RNFRSCal, 1);
};

/* As a sanity check, we put setup parameters in the container, 
 * not in the processor,
 * to have it as an object in the ROOT file for later. */

struct TPCParam {
	constexpr static std::size_t N_PARAMS = 4;

	std::array<double, 2> x_offset;
	std::array<double, 2> x_factor;
	std::array<double, 4> y_offset;
	std::array<double, 4> y_factor;
	
	template<std::size_t N>
	decltype(auto) get() noexcept {
		/* Brackets around return value because of decltype deduction rules. Check link below:
		 * https://stackoverflow.com/questions/27557369/why-does-decltypeauto-return-a-reference-here
		 */
		if constexpr(N == 0)      return (x_offset);
		else if constexpr(N == 1) return (x_offset);
		else if constexpr(N == 2) return (y_offset); 
		else if constexpr(N == 3) return (y_factor); 
		else static_assert(N < 4, "Index out of range.");
	}

	virtual ~TPCParam() = default;
	ClassDef(TPCParam, 1);
};

class TFRSCalCont: public TContainer<RNFRSCal> {
public:
	static nlohmann::json setup; 
	
	TH2I* h1_xy_s2_before_target;
	TH2I* h1_ab_s2_before_target;
	TH2I* h1_xy_s2_after_target;
	TH2I* h1_ab_s2_after_target;

	std::string *setupName{};
	std::array<TPCParam, 7> *tpc_param{};

	TFRSCalCont() ;
	void Init(TDictInfo info) /* override */;
};
