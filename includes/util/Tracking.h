#pragma once

#include <array>
#include "Geometry.h"

class TH2D;

void FillTrack(TH2D* , const std::array<double, 2>& ); 
void FillTrack(TH2D* , const mnd::geom::Line2D& ); 
