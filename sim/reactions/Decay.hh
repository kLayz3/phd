#pragma once

#include "util/MPhysics.h"

#include "G4VDiscreteProcess.hh"
#include "G4VUserTrackInformation.hh"
#include <array>

struct Resonance {
	int id;
	phy::Nucleus mother;
	std::vector<phy::Nucleus> daughters;

	double excitation, width, cumulative_probability;
};

struct DaughterID final : G4VUserTrackInformation {
	explicit DaughterID(size_t i) : index(i) {}
	size_t index;
};

struct ReactionEvent {
	const Resonance* channel = nullptr;
	double requested_z = NAN;
	double mother_beta_gamma = NAN;
	std::vector<std::array<double, 8>> measurements;
};

G4ParticleDefinition* particle(const phy::Nucleus& ion);

class Decay final : public G4VDiscreteProcess {
	ReactionEvent& event;

public:
	explicit Decay(ReactionEvent& e) : G4VDiscreteProcess("PromptDecay"), event(e) {}
	double PostStepGetPhysicalInteractionLength(const G4Track&, double, G4ForceCondition*) override;
	G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
	double GetMeanFreePath(const G4Track&, double, G4ForceCondition*) override { return DBL_MAX; }
};
