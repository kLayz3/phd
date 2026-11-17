#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include "monad/accessories/timing.hxx"

int main(int argc, char** argv) {
    // Let Geant4 choose the MT run manager when the installed Geant4 was
    // built with multithreading support; otherwise this falls back to serial.
	
	srand(time(NULL));
	std::vector<TimePoint> tv;
	tv.emplace_back(TimePoint("start"));

    auto* runManager = G4RunManagerFactory::CreateRunManager();

    auto* detector = new DetectorConstruction;
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(new ActionInitialization(detector));
    const G4String macro = (argc > 1) ? argv[1] : REACTION_DEFAULT_MACRO;
    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + macro);

    delete runManager;

	tv.emplace_back(TimePoint("end"));
	PrintElapsed<kSECOND>(std::move(tv));

    return 0;
}
