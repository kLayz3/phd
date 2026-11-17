#pragma once

#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class C9FragmentationProcess final : public G4VDiscreteProcess {
public:
    C9FragmentationProcess();
    ~C9FragmentationProcess() override = default;

    G4double GetMeanFreePath(
        const G4Track& track,
        G4double previousStepSize,
        G4ForceCondition* condition) override;

    G4VParticleChange* PostStepDoIt(
        const G4Track& track,
        const G4Step& step) override;

private:
    struct Resonance {
        G4int index;
        G4double ex;
        G4double width;
    };

    static G4bool IsC9(const G4Track& track);
    static Resonance SampleResonance(G4double thresholdEx);
    static G4double SampleTruncatedBreitWigner(
        G4double mean, G4double width, G4double thresholdEx);
};
