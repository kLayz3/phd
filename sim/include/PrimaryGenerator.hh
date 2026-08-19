#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include <memory>

class DetectorConstruction;
class G4Event;
class G4ParticleGun;

class PrimaryGenerator final : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGenerator(const DetectorConstruction* detector);
    ~PrimaryGenerator() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    const DetectorConstruction* fDetector;
    std::unique_ptr<G4ParticleGun> fGun;
};
