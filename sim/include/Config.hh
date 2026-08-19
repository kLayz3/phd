#pragma once

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace cfg {

// Beam
inline const G4double beamEnergyPerU = 580.0 * MeV;
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

// IMPORTANT: this is a placeholder TOTAL 9C-removal/reaction cross section.
// It controls the reaction-depth distribution. Replace with the cross section
// appropriate for your beam energy / experiment.
inline const G4double sigmaTotal = 1.0 * barn;

// Effective 8B* model.  These are approximate resonance parameters.
// The relative population is deliberately arbitrary for this first transport test.
inline const G4double ex1 = 0.770 * MeV;
inline const G4double width1 = 35.6 * keV;
inline const G4double ex2 = 2.320 * MeV;
inline const G4double width2 = 350.0 * keV;
inline constexpr G4double state1Probability = 0.5;

// Truncate the Cauchy/Breit-Wigner sampling to avoid very remote tails in this
// simple toy model. This is not an R-matrix line shape.
inline constexpr G4double widthWindow = 5.0;

} // namespace cfg
