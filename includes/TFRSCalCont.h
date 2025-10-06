#pragma once

#include "TContainer.hxx"
#include "nlohmann/json.hpp"
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

class TFRSCalCont: public TContainer<RNFRSCal> {
public:
	using TPCParam = std::tuple<
		std::array<double, 2>, // x_offset
		std::array<double, 2>, // x_factor
		std::array<double, 4>, // y_offset
		std::array<double, 4>  // y_factor
	>;

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
