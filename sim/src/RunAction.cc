#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction() {
    auto* a = G4AnalysisManager::Instance();
    a->SetVerboseLevel(1);
    a->SetDefaultFileType("root");
    a->SetFileName("c9frag");
    a->SetNtupleMerging(true);

    a->CreateNtuple("events", "9C -> 8B* -> p + 7Be transport");
    a->CreateNtupleIColumn("event");              // 0
    a->CreateNtupleIColumn("reacted");            // 1
    a->CreateNtupleIColumn("resonance");          // 2
    a->CreateNtupleDColumn("reaction_z_mm");      // 3
    a->CreateNtupleDColumn("depth_mg_cm2");       // 4
    a->CreateNtupleDColumn("c9_T_per_u_MeV");     // 5
    a->CreateNtupleDColumn("Ex_MeV");             // 6
    a->CreateNtupleDColumn("theta_vertex_deg");   // 7
    a->CreateNtupleDColumn("theta_exit_deg");     // 8
    a->CreateNtupleIColumn("both_exited");        // 9
    a->CreateNtupleDColumn("p_vertex_MeVc");      // 10
    a->CreateNtupleDColumn("be7_vertex_MeVc");    // 11
    a->CreateNtupleDColumn("p_exit_MeVc");        // 12
    a->CreateNtupleDColumn("be7_exit_MeVc");      // 13
    a->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
    G4AnalysisManager::Instance()->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto* a = G4AnalysisManager::Instance();
    a->Write();
    a->CloseFile();
}
