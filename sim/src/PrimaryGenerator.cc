#include "PrimaryGenerator.hh"

#include "Config.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

PrimaryGenerator::PrimaryGenerator(const DetectorConstruction* detector)
    : fDetector(detector), fGun(std::make_unique<G4ParticleGun>(1)) {}

PrimaryGenerator::~PrimaryGenerator() = default;

void PrimaryGenerator::GeneratePrimaries(G4Event* event) {
    auto* ion = G4IonTable::GetIonTable()->GetIon(cfg::beamZ, cfg::beamA, 0.0);

    fGun->SetParticleDefinition(ion);
    fGun->SetParticleCharge(cfg::beamZ * eplus); // 9C6+
    fGun->SetParticleEnergy(cfg::beamA * cfg::beamEnergyPerU);
    fGun->SetParticleMomentumDirection({0.0, 0.0, 1.0});
    fGun->SetParticlePosition({0.0, 0.0, fDetector->TargetEntranceZ() - 1.0 * um});
    fGun->GeneratePrimaryVertex(event);
}
