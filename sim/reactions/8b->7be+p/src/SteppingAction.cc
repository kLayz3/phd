#include "SteppingAction.hh"

#include "EventInfo.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    const auto* prePV = pre->GetPhysicalVolume();
    const auto* postPV = post->GetPhysicalVolume();

    if (prePV == nullptr || postPV == nullptr)
        return;
    if (prePV->GetName() != "TargetPV" || postPV->GetName() != "WorldPV")
        return;
    if (post->GetStepStatus() != fGeomBoundary)
        return;

    const auto* track = step->GetTrack();
    const auto* creator = track->GetCreatorProcess();
    if (creator == nullptr || creator->GetProcessName() != "C9Fragmentation")
        return;

    const auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    if (event == nullptr)
        return;
    auto* info = static_cast<EventInfo*>(event->GetUserInformation());
    if (info == nullptr)
        return;

    const auto* def = track->GetParticleDefinition();
    const G4ThreeVector momentum = post->GetMomentum();

    if (def == G4Proton::ProtonDefinition()) {
        info->protonExited = true;
        info->protonExitMomentum = momentum;
        return;
    }

    auto* be7 = G4IonTable::GetIonTable()->GetIon(4, 7, 0.0);
    if (def == be7) {
        info->be7Exited = true;
        info->be7ExitMomentum = momentum;
    }
}
