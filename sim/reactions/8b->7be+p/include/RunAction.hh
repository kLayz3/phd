#pragma once

#include "G4UserRunAction.hh"

class G4Run;

class RunAction final : public G4UserRunAction {
public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
};
