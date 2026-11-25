#pragma once
#include "TransportConfig.hh"
class G4LogicalVolume;
namespace transport {
void place_telescope(const Config&, G4LogicalVolume* world, double half_width, double origin_z, bool use_vacuum);
}
