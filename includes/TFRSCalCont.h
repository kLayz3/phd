#pragma once

#include "TContainer.h"
#include "nlohmann/json.hpp"

class TFRSCalCont: public TContainer {
public:
	class ToF {
	public:
		double val; // ns;
		ClassDef(ToF, 1);
	};

	class Position {
	public:
		double x{}; double y{};
		double a{}; double b{};
		Position() {};
		virtual ~Position() = default;
		ClassDef(Position, 1);
	};
	
	static nlohmann::json setup; 
	Position before_target{}, after_target{};

	TFRSCalCont();
	~TFRSCalCont();
	ClassDef(TFRSCalCont, 1);
};
