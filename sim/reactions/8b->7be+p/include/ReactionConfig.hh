#pragma once

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace reaction {

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
inline constexpr G4double state1Probability = 0.18;

// Truncate the Cauchy/Breit-Wigner sampling to avoid very remote tails in this
// simple toy model. This is not an R-matrix line shape.
inline constexpr G4double widthWindow = 5.0;

} // namespace reaction
