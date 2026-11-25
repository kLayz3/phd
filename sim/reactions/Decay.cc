#include "Decay.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TGenPhaseSpace.h"
#include <numeric>

G4ParticleDefinition* particle(const phy::Nucleus& ion) {
	if(ion.A == 1 && ion.Z == 1) return G4Proton::Definition();
	
	auto* def = G4IonTable::GetIonTable()->GetIon(ion.Z, ion.A, 0.);
	if(!def)
		ERROR("Geant4 has no ion with Z={}, A={}", ion.Z, ion.A);

	return def;
}

double Decay::PostStepGetPhysicalInteractionLength(const G4Track& t, double, G4ForceCondition* condition) {
	*condition = NotForced;
	if(!event.channel || t.GetParentID() || t.GetVolume()->GetName() != "Target")
		return DBL_MAX;
	const double dz = event.requested_z - t.GetPosition().z();
	const double direction = t.GetMomentumDirection().z();
	
	if(direction <= 0)
		return DBL_MAX;
	
	return std::max(0., dz / direction);
}

G4VParticleChange* Decay::PostStepDoIt(const G4Track& t, const G4Step&) {
	aParticleChange.Initialize(t);
	const auto& channel = *event.channel;
	const double ground = particle(channel.mother)->GetPDGMass();
	std::vector<double> masses;
	for(const auto& ion: channel.daughters) masses.push_back(particle(ion)->GetPDGMass());
	const double threshold = std::accumulate(masses.begin(), masses.end(), 0.) - ground;
	double ex = channel.excitation;
	if(channel.width > 0) {
		const double half = channel.width / 2;
		const double low = std::max(threshold, ex - 5 * channel.width);
		const double high = ex + 5 * channel.width;
		if(low >= high) ERROR("Reaction {} has no open decay window", channel.id);
		const double a = atan((low-ex)/half), b = atan((high-ex)/half);
		ex += half * tan(a + (b-a)*G4UniformRand());
	}
	if(ex <= threshold) ERROR("Reaction {} is below threshold", channel.id);
	const double kinetic = t.GetKineticEnergy();
	const double mass = ground + ex, energy = kinetic + mass;
	const auto p = sqrt(kinetic * (kinetic + 2*mass)) * t.GetMomentumDirection();
	event.mother_beta_gamma = p.mag() / mass;
	TLorentzVector mother(p.x(), p.y(), p.z(), energy);
	TGenPhaseSpace decay;
	if(!decay.SetDecay(mother, masses.size(), masses.data())) ERROR("Invalid decay {}", channel.id);
	// Generate uniform Lorentz-invariant phase space (Generate returns a normalized weight).
	while(decay.Generate() < G4UniformRand()) {}
	aParticleChange.SetNumberOfSecondaries(masses.size());
	for(size_t i = 0; i < masses.size(); ++i) {
		const auto& v = *decay.GetDecay(i);
		auto* dp = new G4DynamicParticle(particle(channel.daughters[i]), G4ThreeVector(v.Px(),v.Py(),v.Pz()));
		auto* secondary = new G4Track(dp, t.GetGlobalTime(), t.GetPosition());
		secondary->SetUserInformation(new DaughterID(i));
		aParticleChange.AddSecondary(secondary);
	}
	aParticleChange.ProposeTrackStatus(fStopAndKill);
	aParticleChange.ProposeEnergy(0.);
	return &aParticleChange;
}
