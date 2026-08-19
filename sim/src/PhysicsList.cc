#include "PhysicsList.hh"

#include "FragmentationPhysics.hh"
#include "Config.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() {
    SetVerboseLevel(1);
    SetDefaultCutValue(cfg::defaultProductionCut);

    // EM only: ionisation, energy-loss fluctuations, multiple scattering, etc.
    RegisterPhysics(new G4EmStandardPhysics_option4(1));

    // Exactly one nuclear reaction model, controlled by us.
    RegisterPhysics(new FragmentationPhysics());
}

void PhysicsList::SetCuts() {
    // Establish the default production cuts first, then override electrons.
    // This is a production threshold, not a tracking cut: energy transfers
    // below the threshold remain part of continuous dE/dx.
    SetCutsWithDefault();
    SetCutValue(cfg::deltaElectronRangeCut, "e-");
}
