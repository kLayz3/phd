#pragma once

#include "G4UserEventAction.hh"

class DetectorConstruction;
class G4Event;

class EventAction final : public G4UserEventAction {
public:
    explicit EventAction(const DetectorConstruction* detector);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

private:
    const DetectorConstruction* fDetector;
};
