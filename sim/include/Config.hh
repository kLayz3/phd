#pragma once

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace cfg {

// EM production cuts.  Keep the usual fine cut for everything except
// secondary electrons: a large e- range cut strongly suppresses explicit
// delta-ray tracks while retaining their sub-threshold energy transfer in
// the continuous ionisation energy loss.
inline const G4double defaultProductionCut = 0.1 * mm;
inline const G4double deltaElectronRangeCut = 10.0 * mm;

// Shared world/target dimensions and beam entrance state.
inline const G4double worldHalfXY = 25.0 * cm;
inline const G4double beamStartGap = 1.0 * um;

} // namespace cfg
