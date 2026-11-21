#include "../../includes/monad/accessories/formatting.hxx"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TVirtualPad.h>
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

constexpr auto nil = std::nullopt;

// root -l 'sim/transport/draw_transport.C("transport_9C.root")'
void draw_transport(const char* filename = "foot_transport.root", const char* pdf = "") {
	std::unique_ptr<TFile> file{TFile::Open(filename, "READ")};
	if(!file || file->IsZombie()) {
		ERROR("Cannot open {}", filename);
	}
	auto* tree = file->Get<TTree>("crossings");
	if(!tree || tree->GetEntries() == 0) {
		ERROR("Missing or empty crossings tree in {}", filename);
	}
	static u32 counter = 0;
	auto draw = [&](const std::string& expression, const std::string& cut,
		const std::string& title,
		std::optional<u32> nbinsx = nil,
		std::optional<std::array<double,2>> xrange = nil,
		bool two_dimensional = false,
		std::optional<u32> nbinsy = nil,
		std::optional<std::array<double,2>> yrange = nil)
	{
		const std::string name = "transport_plot_" + std::to_string(counter++);
		std::string draw_expr = expression + ">>" + name;
		draw_expr += "(";
		if(nbinsx.has_value()) {
			draw_expr += fmt::to_string(*nbinsx);
		}
		if(xrange.has_value()) {
			if(!nbinsx.has_value())
				ERROR("Must give also nbinsx, hello?");
			draw_expr += fmt::format(",{:.6f},{:.6f}", (*xrange)[0], (*xrange)[1]);
		}
		if(two_dimensional) {
			if(nbinsy.has_value()) {
				if(!nbinsx.has_value() or !xrange.has_value()) /* In this case, x ranges must also be given! */
					ERROR("Wanting to pass in nbinsy, but nbinsx/xrange not fully quantified.");
				draw_expr += fmt::format(",{}", *nbinsy);
			}
			if(yrange.has_value()) {
				if(!nbinsy.has_value())
					ERROR("Must give also nbinsy, hello?");
				draw_expr += fmt::format(",{:.6f},{:.6f}", (*xrange)[0], (*xrange)[1]);
			}
		}
		draw_expr += ")";
		if(tree->Draw(draw_expr.c_str(), cut.c_str(), "goff") < 0) {
			ERROR("Cannot plot {}", expression);
		}
		std::unique_ptr<TH1> h{tree->GetHistogram()};
		if(!h) ERROR("Cannot retrieve histogram for {}", expression);
		h->SetDirectory(nullptr);
		h->SetTitle(title.c_str());
		h->SetFillColor(kRed+2);
		h->SetFillStyle(1001);
		// Pads retain their own histograms after this helper and the input file are destroyed.
		h->DrawCopy(two_dimensional ? "COLZ" : "HIST");
		gPad->SetGrid();
	};

	auto* measurements = new TCanvas("TransportMeasurements", "Unsmeared scalar measurements", 1800, 900);
	measurements->Divide(4, 2);

	auto* positions = new TCanvas("TransportPositions", "True midplane positions", 1800, 900);
	positions->Divide(4, 2);

	auto* angles = new TCanvas("TransportDeflections", "Accumulated deflections from entrance direction", 1800, 1400);
	angles->Divide(4, 4);

	constexpr u32 nbins = 300;
	std::string number_present;
	for(u32 i=0; i<8; ++i) {
		const std::string prefix = "d" + std::to_string(i) + "_";
		const std::string label = "Detector " + std::to_string(i);
		const std::string cut = prefix + "present";
		if(i) number_present += "+";
		number_present += cut;
		measurements->cd(i+1);
		draw(prefix+"measurement_mm", cut, label+";Bare readout-axis position [mm];Count");

		positions->cd(i+1);
		draw(prefix+"y_mm:"+prefix+"x_mm", cut, label+";True x [mm];True y [mm]",
			nbins, {{-3,3}}, /*2D*/true,
			nbins, {{-3,3}}
		);

		angles->cd(i+1);
		draw("1000*("+prefix+"ax-beam_ax)", cut, label+";a_{x} - a_{x,entrance} [mrad];Count",
			nbins, {{-4, 4}} );
		angles->cd(i+9);
		draw("1000*("+prefix+"ay-beam_ay)", cut, label+";a_{y} - a_{y,entrance} [mrad];Count",
			nbins, {{-4, 4}} );
	}

	auto* overview = new TCanvas("TransportOverview", "Telescope transport overview", 1400, 1000);
	overview->Divide(2, 2);
	overview->cd(1);
	draw(number_present, "", "Crossed measurement planes;Number of planes;Events");
	overview->cd(2);
	draw("1000*(d7_ax-beam_ax)", "d7_present", "Last-plane deflection;a_{x} - a_{x,entrance} [mrad];Count",
	  nbins, {{-3,3}});
	overview->cd(3);
	draw("1000*(d7_ay-beam_ay):1000*(d7_ax-beam_ax)", "d7_present",
		"Last-plane deflection;#Delta a_{x} [mrad];#Delta a_{y} [mrad]",
			nbins, {{-3,3}}, /*2D*/true,
			nbins, {{-3,3}}
	);

	overview->cd(4);
	draw("d0_kinetic_MeV-d7_kinetic_MeV", "d0_present && d7_present",
		"Energy lost between first and last midplanes;Energy loss [MeV];Count");

	const std::vector<TCanvas*> canvases{overview, measurements, positions, angles};
	for(auto* c: canvases) {
		c->Modified();
		c->Update();
	}

	if(pdf && *pdf) {
		WARN("Dumping the canvases into a pdf: '{}'", pdf);
		overview->Print((std::string(pdf)+"[").c_str());
		for(auto* c: canvases) c->Print(pdf);
		overview->Print((std::string(pdf)+"]").c_str());
	}

	WARN("{} generated particles. Positions include transport but no detector smearing.\n", tree->GetEntries());
}
