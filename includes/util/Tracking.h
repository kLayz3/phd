#pragma once

#include <array>
#include "Geometry.h"

class TH2D;

std::array<double, 2> GetLine(const mnd::geom::Point& , const mnd::geom::Point& ) noexcept;
std::array<double, 2> GetLine(const std::array<std::array<double, 2>, 2> ) noexcept;
void FillTrack(TH2D* , const std::array<double, 2>& a); 
