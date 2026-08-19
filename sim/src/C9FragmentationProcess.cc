#include "C9FragmentationProcess.hh"

#include "Config.hh"
#include "EventInfo.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <limits>

C9FragmentationProcess::C9FragmentationProcess()
    : G4VDiscreteProcess("C9Fragmentation") {}

G4bool C9FragmentationProcess::IsC9(const G4Track& track) {
    const auto* def = track.GetParticleDefinition();
    return def->GetAtomicNumber() == cfg::beamZ
        && def->GetAtomicMass() == cfg::beamA;
}

G4double C9FragmentationProcess::GetMeanFreePath(
    const G4Track& track,
    G4double,
    G4ForceCondition* condition)
{
    *condition = NotForced;

    if (!IsC9(track))
        return std::numeric_limits<G4double>::max();

    const auto* material = track.GetMaterial();
    if (material == nullptr || material->GetName() != "Be9Target")
        return std::numeric_limits<G4double>::max();

    // For this pure-element target, n * sigma is the macroscopic cross section.
    const G4double nAtoms = material->GetTotNbOfAtomsPerVolume();
    return 1.0 / (nAtoms * cfg::sigmaTotal);
}

G4double C9FragmentationProcess::SampleTruncatedBreitWigner(
    G4double mean,
    G4double width,
    G4double thresholdEx)
{
    // Non-relativistic Breit-Wigner == Cauchy distribution with scale Gamma/2.
    // The finite window is deliberate: this is a transport toy model, not an
    // R-matrix treatment of the resonance line shape.
    const G4double lo = std::max(
        thresholdEx, mean - cfg::widthWindow * width);
    const G4double hi = mean + cfg::widthWindow * width;

    while (true) {
        const G4double u = G4UniformRand();
        const G4double x = mean + 0.5 * width * std::tan(pi * (u - 0.5));
        if (x >= lo && x <= hi)
            return x;
    }
}

C9FragmentationProcess::Resonance
C9FragmentationProcess::SampleResonance(G4double thresholdEx) {
    if (G4UniformRand() < cfg::state1Probability) {
        return {1, SampleTruncatedBreitWigner(
            cfg::ex1, cfg::width1, thresholdEx)};
    }
    return {2, SampleTruncatedBreitWigner(
        cfg::ex2, cfg::width2, thresholdEx)};
}

G4VParticleChange* C9FragmentationProcess::PostStepDoIt(
    const G4Track& track,
    const G4Step&)
{
    aParticleChange.Initialize(track);

    auto* ionTable = G4IonTable::GetIonTable();
    auto* b8gs = ionTable->GetIon(5, 8, 0.0);
    auto* be7 = ionTable->GetIon(4, 7, 0.0);
    auto* proton = G4Proton::ProtonDefinition();

    const G4double mB8gs = b8gs->GetPDGMass();
    const G4double mBe7 = be7->GetPDGMass();
    const G4double mp = proton->GetPDGMass();
    const G4double thresholdEx = mBe7 + mp - mB8gs;

    const Resonance resonance = SampleResonance(thresholdEx);
    const G4double parentMass = mB8gs + resonance.ex;

    // Effective projectile-fragment production model:
    // 8B* keeps the instantaneous 9C velocity at the reaction vertex.
    // The unobserved target/removal system X absorbs the missing 4-momentum.
    const G4ThreeVector beta = track.GetMomentum() / track.GetTotalEnergy();
    const G4double beta2 = beta.mag2();
    const G4double gamma = 1.0 / std::sqrt(1.0 - beta2);
    G4LorentzVector pB8(gamma * parentMass * beta, gamma * parentMass);

    // Exact relativistic two-body momentum in the 8B* rest frame.
    const G4double M2 = parentMass * parentMass;
    const G4double sum2 = (mp + mBe7) * (mp + mBe7);
    const G4double diff2 = (mp - mBe7) * (mp - mBe7);
    const G4double lambda = (M2 - sum2) * (M2 - diff2);
    const G4double q = std::sqrt(std::max(0.0, lambda)) / (2.0 * parentMass);

    // Isotropic decay in the 8B* frame. Alignment can be inserted here later.
    const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
    const G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    const G4double phi = twopi * G4UniformRand();
    const G4ThreeVector nhat(
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta);

    G4LorentzVector pPStar(q * nhat, std::sqrt(mp * mp + q * q));
    G4LorentzVector pBeStar(-q * nhat, std::sqrt(mBe7 * mBe7 + q * q));

    const G4ThreeVector boost = pB8.vect() / pB8.e();
    pPStar.boost(boost);
    pBeStar.boost(boost);

    auto* pDynamic = new G4DynamicParticle(
        proton, pPStar.vect().unit(), pPStar.e() - mp);
    auto* beDynamic = new G4DynamicParticle(
        be7, pBeStar.vect().unit(), pBeStar.e() - mBe7);

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(pDynamic);
    aParticleChange.AddSecondary(beDynamic);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.ProposeEnergy(0.0);
    aParticleChange.ProposeLocalEnergyDeposit(0.0);

    // Save generator-level truth for this event.
    const auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    if (event != nullptr) {
        auto* info = static_cast<EventInfo*>(event->GetUserInformation());
        if (info != nullptr) {
            info->reacted = true;
            info->resonance = resonance.index;
            info->reactionZ = track.GetPosition().z();
            info->reactionTrackLength = track.GetTrackLength();
            info->c9KineticPerU = track.GetKineticEnergy() / cfg::beamA;
            info->excitationEnergy = resonance.ex;
            info->protonVertexMomentum = pPStar.vect();
            info->be7VertexMomentum = pBeStar.vect();
        }
    }

    return &aParticleChange;
}
