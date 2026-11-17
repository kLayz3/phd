#pragma once

#include "G4ThreeVector.hh"
#include "G4VUserEventInformation.hh"
#include "globals.hh"

class EventInfo final : public G4VUserEventInformation {
public:
    void Print() const override {}

    G4bool reacted = false;
    G4int resonance = -1;

    G4double reactionZ = 0.0;
    G4double reactionTrackLength = 0.0;
    G4double c9KineticPerU = 0.0;
    G4double excitationEnergy = 0.0;

    G4ThreeVector protonVertexMomentum;
    G4ThreeVector be7VertexMomentum;

    G4bool protonExited = false;
    G4bool be7Exited = false;
    G4ThreeVector protonExitMomentum;
    G4ThreeVector be7ExitMomentum;
};
