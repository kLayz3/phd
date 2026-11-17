#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "EventInfo.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <limits>

EventAction::EventAction(const DetectorConstruction* detector)
    : fDetector(detector) {}

void EventAction::BeginOfEventAction(const G4Event*) {
    G4EventManager::GetEventManager()->SetUserInformation(new EventInfo);
}

void EventAction::EndOfEventAction(const G4Event* event) {
    const auto* info = static_cast<const EventInfo*>(event->GetUserInformation());
    auto* a = G4AnalysisManager::Instance();

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const bool reacted = info != nullptr && info->reacted;
    const bool bothExited = reacted && info->protonExited && info->be7Exited;

    double reactionZ = nan;
    double depth = nan;
    double tPerU = nan;
    double ex = nan;
    double thetaVertex = nan;
    double thetaExit = nan;
    double pVertex = nan;
    double beVertex = nan;
    double pExit = nan;
    double beExit = nan;
    int resonance = -1;

    if (reacted) {
        reactionZ = info->reactionZ / mm;
        const G4double longitudinalDepth =
            info->reactionZ - fDetector->TargetEntranceZ();
        const G4double arealDepth = longitudinalDepth * fDetector->TargetDensity();
        depth = arealDepth / (mg / (cm * cm));
        tPerU = info->c9KineticPerU / MeV;
        ex = info->excitationEnergy / MeV;
        thetaVertex = info->protonVertexMomentum.angle(info->be7VertexMomentum) / deg;
        pVertex = info->protonVertexMomentum.mag() / MeV;
        beVertex = info->be7VertexMomentum.mag() / MeV;

        if (bothExited) {
            thetaExit = info->protonExitMomentum.angle(info->be7ExitMomentum) / deg;
            pExit = info->protonExitMomentum.mag() / MeV;
            beExit = info->be7ExitMomentum.mag() / MeV;
        }
        resonance = info->resonance;
    }

    a->FillNtupleIColumn(0, event->GetEventID());
    a->FillNtupleIColumn(1, reacted ? 1 : 0);
    a->FillNtupleIColumn(2, resonance);
    a->FillNtupleDColumn(3, reactionZ);
    a->FillNtupleDColumn(4, depth);
    a->FillNtupleDColumn(5, tPerU);
    a->FillNtupleDColumn(6, ex);
    a->FillNtupleDColumn(7, thetaVertex);
    a->FillNtupleDColumn(8, thetaExit);
    a->FillNtupleIColumn(9, bothExited ? 1 : 0);
    a->FillNtupleDColumn(10, pVertex);
    a->FillNtupleDColumn(11, beVertex);
    a->FillNtupleDColumn(12, pExit);
    a->FillNtupleDColumn(13, beExit);
    a->AddNtupleRow();
}
