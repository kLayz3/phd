#pragma once

#include "G4VPhysicsConstructor.hh"

class FragmentationPhysics final : public G4VPhysicsConstructor {
public:
    FragmentationPhysics();
    ~FragmentationPhysics() override = default;

    void ConstructParticle() override;
    void ConstructProcess() override;
};
