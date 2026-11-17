#include "FragmentationPhysics.hh"

#include "C9FragmentationProcess.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"

FragmentationPhysics::FragmentationPhysics()
    : G4VPhysicsConstructor("FragmentationPhysics") {}

void FragmentationPhysics::ConstructParticle() {
    G4GenericIon::GenericIonDefinition();
    G4Proton::ProtonDefinition();
}

void FragmentationPhysics::ConstructProcess() {
    // This is a user-defined G4VDiscreteProcess. Registering it through
    // G4PhysicsListHelper would require a process type/subtype present in
    // the helper's ordering table.  For a simple PostStep-only process,
    // registering directly with the particle's process manager is the
    // appropriate Geant4 API.
    auto* genericIon = G4GenericIon::GenericIonDefinition();
    auto* processManager = genericIon->GetProcessManager();
    processManager->AddDiscreteProcess(new C9FragmentationProcess);
}
