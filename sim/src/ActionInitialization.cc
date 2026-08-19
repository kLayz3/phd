#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* detector)
    : fDetector(detector) {}

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction);
}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGenerator(fDetector));
    SetUserAction(new RunAction);
    SetUserAction(new EventAction(fDetector));
    SetUserAction(new SteppingAction);
}
