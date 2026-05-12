#pragma once

#include "PolyFitter.hxx"


/* One extra algorithm to solve general 2D linear problem,
 * arising from solving rotational measurements:
 * x' = cos(t)*x + sin(t)*y  , AKA:
 * x' = a*x + b*y      where a^2 + b^2 == 1
 * Where `(x,y)` is the 'true' referent measurement,
 * and `x'` is what the detector gives us.
 * We want to solve for `a` and `b`.
 * We are given sequences of events: `(x,y, x')` here
 * given as the vectors `x0`, `y0`, `x` respectively.
 */

struct AngleFitResult {
	enum class Direction { X, Y };
	double a; // cos(tx)
	double b; // sin(tx)
	double Angle(const Direction ) const noexcept; // tx
};

AngleFitResult FitAngle (
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
);

struct AngleOffsetFitResult {
	AngleFitResult t;
	double c; // -dx*cos(t.tx) - dy*sin(t.tx)
};

AngleOffsetFitResult FitAngleOffset (
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
);
