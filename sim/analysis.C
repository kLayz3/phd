#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TROOT.h>

void analysis(const char* file = "c9frag.root") {
    ROOT::RDataFrame df("events", file);
    auto d = df.Filter("reacted && both_exited");

    auto hDepth = d.Histo1D(
        {"hDepth", ";reaction depth [mg/cm^{2}];events", 100, 0, 7500},
        "depth_mg_cm2");

    auto hCorr = d.Histo2D(
        {"hCorr", ";#theta_{vertex} [deg];#theta_{exit} [deg]",
         160, 0, 8, 160, 0, 8},
        "theta_vertex_deg", "theta_exit_deg");

    auto hDepthAngle = d.Histo2D(
        {"hDepthAngle", ";reaction depth [mg/cm^{2}];#theta_{exit} [deg]",
         100, 0, 7500, 160, 0, 8},
        "depth_mg_cm2", "theta_exit_deg");

    auto h1 = d.Filter("resonance == 1").Histo1D(
        {"h1", "0.77 MeV state;#theta_{exit} [deg];events", 160, 0, 8},
        "theta_exit_deg");
    auto h2 = d.Filter("resonance == 2").Histo1D(
        {"h2", "2.32 MeV state;#theta_{exit} [deg];events", 160, 0, 8},
        "theta_exit_deg");

    auto* c1 = new TCanvas("c_depth", "reaction depth", 900, 700);
    hDepth->Draw();

    auto* c2 = new TCanvas("c_corr", "vertex vs exit", 900, 700);
    hCorr->Draw("COLZ");

    auto* c3 = new TCanvas("c_depth_angle", "depth vs angle", 900, 700);
    hDepthAngle->Draw("COLZ");

    auto* c4 = new TCanvas("c_states", "states", 900, 700);
    h1->SetLineColor(kBlue + 1);
    h2->SetLineColor(kRed + 1);
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
}
