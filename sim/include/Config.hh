#pragma once

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace cfg {

// Beam
//inline const G4double beamEnergyPerU = 580.0 * MeV;
inline const G4double beamEnergyPerU = 700.0 * MeV;
inline const G4int beamA = 9;
inline const G4int beamZ = 6;

// Target: 7500 mg/cm^2 of Be, using 1.848 g/cm^3.
inline const G4double targetArealDensity = 7500.0 * mg / (cm * cm);
inline const G4double targetDensity = 1.848 * g / (cm * cm * cm);

// EM production cuts.  Keep the usual fine cut for everything except
// secondary electrons: a large e- range cut strongly suppresses explicit
// delta-ray tracks while retaining their sub-threshold energy transfer in
// the continuous ionisation energy loss.
inline const G4double defaultProductionCut = 0.1 * mm;
inline const G4double deltaElectronRangeCut = 10.0 * mm;

// Shared world/target dimensions and beam entrance state.
inline const G4double worldHalfXY = 25.0 * cm;
inline const G4double worldHalfZ = 50.0 * cm;
inline const G4double targetHalfXY = 10.0 * cm;
inline const G4double beamX = 0.0 * mm;
inline const G4double beamY = 0.0 * mm;
inline const G4double beamAx = 0.0;
inline const G4double beamAy = 0.0;
inline const G4double beamStartGap = 1.0 * um;

} // namespace cfg
