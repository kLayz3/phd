#pragma once

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization final : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(const DetectorConstruction* detector);

    void BuildForMaster() const override;
    void Build() const override;

private:
    const DetectorConstruction* fDetector;
};
