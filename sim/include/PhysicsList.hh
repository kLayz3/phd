#pragma once

#include "G4VModularPhysicsList.hh"

class PhysicsList final : public G4VModularPhysicsList {
public:
    PhysicsList();
    ~PhysicsList() override = default;

    void SetCuts() override;
};
