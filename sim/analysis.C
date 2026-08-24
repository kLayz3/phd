#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TROOT.h>

constexpr double THETA_VERTEX_MAX = 100;
constexpr double THETA_EXIT_MAX   = 100;

void analysis(const char* file = "c9frag.root") {
    ROOT::RDataFrame df("events", file);
    auto d = df
		.Filter("reacted && both_exited")
		.Define("theta_vertex_mrad", "theta_vertex_deg * M_PI * 1000.0/180")
		.Define("theta_exit_mrad",   "theta_exit_deg   * M_PI * 1000.0/180");

    auto hDepth = d.Histo1D(
        {"hDepth", ";reaction depth [mg/cm^{2}];events", 100, 0, 7500},
        "depth_mg_cm2");

    auto hCorr = d.Histo2D(
        {"hCorr", ";#theta_{vertex} [mrad];#theta_{exit} [mrad]",
         160, 0, THETA_VERTEX_MAX,
		 160, 0, THETA_EXIT_MAX},
        "theta_vertex_mrad", "theta_exit_mrad");

    auto hDepthAngle = d.Histo2D(
        {"hDepthAngle", ";reaction depth [mg/cm^{2}];#theta_{exit} [mrad]",
         100, 0, 7500, 160, 0, THETA_EXIT_MAX},
        "depth_mg_cm2", "theta_exit_mrad");

    auto h1 = d.Filter("resonance == 1").Histo1D(
        {"h1", "0.77 MeV state;#theta_{exit} [mrad];events", 160, 0, THETA_EXIT_MAX},
        "theta_exit_mrad");
    auto h2 = d.Filter("resonance == 2").Histo1D(
        {"h2", "2.32 MeV state;#theta_{exit} [mrad];events", 160, 0, THETA_EXIT_MAX},
        "theta_exit_mrad");
    auto hSum = d.Histo1D(
        {"hSum", "both states;#theta_{exit} [mrad];events", 160, 0, THETA_EXIT_MAX},
        "theta_exit_mrad");

    auto* c1 = new TCanvas("c_depth", "reaction depth", 900, 700);
    hDepth->DrawCopy();

    auto* c2 = new TCanvas("c_corr", "vertex vs exit", 900, 700);
    hCorr->DrawCopy("COLZ");

    auto* c3 = new TCanvas("c_depth_angle", "depth vs angle", 900, 700);
    hDepthAngle->DrawCopy("COLZ");

    auto* c4 = new TCanvas("c_states", "states", 1300, 1300);
	c4->Divide(1,2);
    h1->SetLineColor(kBlue + 1);
	h2->SetLineColor(kRed + 1);
	const double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
	h1->SetMaximum(1.10 * ymax);

	c4->cd(1);
    h1->DrawCopy("HIST");
    h2->DrawCopy("HIST SAME");

	c4->cd(2);
	hSum->DrawCopy("HIST");
}
