#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    G4VPhysicalVolume* Construct() override;

    G4double TargetThickness() const { return fTargetThickness; }
    G4double TargetDensity() const { return fTargetDensity; }
    G4double TargetEntranceZ() const { return -0.5 * fTargetThickness; }

private:
    G4double fTargetThickness;
    G4double fTargetDensity;
};
