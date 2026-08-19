#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

int main(int argc, char** argv) {
    // Let Geant4 choose the MT run manager when the installed Geant4 was
    // built with multithreading support; otherwise this falls back to serial.
    auto* runManager = G4RunManagerFactory::CreateRunManager();

    auto* detector = new DetectorConstruction;
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(new ActionInitialization(detector));
    const G4String macro = (argc > 1) ? argv[1] : "run.mac";
    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + macro);

    delete runManager;
    return 0;
}
