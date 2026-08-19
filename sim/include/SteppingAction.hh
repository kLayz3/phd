#pragma once

#include "G4UserSteppingAction.hh"

class G4Step;

class SteppingAction final : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override;
};
