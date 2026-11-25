#include "Analysis.hh"
#include "Eigen/QR"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include "util/TextUtils.h"

using mnd::geom::Line3D;
using namespace mnd::col::literals;
using namespace std::literals;

namespace {
constexpr auto n_pairs = transport::N / 2;

size_t common_track_count(const std::vector<Resonance>& channels) {
	const size_t count = channels.front().daughters.size();
	return std::all_of(channels.begin(), channels.end(), [count](const auto& r) {
		return r.daughters.size() == count;
	}) ? count : 0;
}
}

void HistConfig::add_options(CLI::App& app) {
	add_logged_option(app, "-o,--save", save, "Save the resulting canvases as one or more extensions.")
		->delimiter(',');
	add_logged_option(app, "-i,--info", info,
		"Sequence of strings, such that the output will be saved descending these directories in autosave.")
		->delimiter(',');
	add_logged_option(app, "--rho-binning", rho_binning, "Binning of the ρ-value, in mrad.")
		->delimiter(',');
	add_logged_flag(app, "--scaling", do_scaling,
		"Rescale ρ using the simulated mother's βγ at the decay point.");
	add_logged_option(app, "--q-binning", q_binning, "Binning of the FOOT tracks' charge values.")
		->delimiter(',');
	add_logged_option(app, "--total-q-binning", all_q_binning,
		"Binning of the wholesome FOOT tracks' charge values. "
		"Done before any FOOT selection; events just need to pass SCI selection.")
		->delimiter(',');
	add_logged_option(app, "--theta-binning", theta_binning,
		"Binning of the θ values, in mrad. This represents the angle between heavy fragment and all the lighter ones.")
		->delimiter(',');
	add_logged_option(app, "--rho-tstyle", rho_style, "Configure the ROOT style of ρ-value 1D histogram.\n"
		"Arg[0]: font style.\nArg[1]: font width (in inches).\nArg[2]: true => show statistics box.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "-d,--distance", distance_cut,
		"Distance from vertex cut (all tracks must be below this threshold). Default no cut.")
		->check(CLI::Range(0.0, 10.0));
}

Histos::ThetaHistograms::ThetaHistograms(const HistConfig& cfg)
	: theta1_vs_theta2{"#theta(p_{1}, HI) [mrad]:#theta(p_{2}, HI) [mrad]@Traced by the FOOT",
		cfg.theta_binning, cfg.theta_binning},
	  rho_vs_opening{"#rho [mrad]:#theta(p_{1}, p_{2}) [mrad]@Traced by the FOOT",
		cfg.theta_binning, cfg.rho_binning} {}

Histos::Histos(const HistConfig& cfg, size_t track_count)
	:
	h1_track_mult{"Track multiplicity [unique tracks]@-1 means no FOOT in event.", kRed-1, 11, -1.5, 9.5},
	h1_foots_fired{"Number of detector pairs fired@0..=4", kBlue-1,
		n_pairs+1, -0.5, n_pairs+0.5},
	h2_q_vs_mult{"Track charge [Q]:Track multiplicity@Full FOOT system",
		10, -0.5, 9.5, 40, 0,8},
	h2_score_vs_mult{"Track score [a.u.]:Track multiplicity@Full FOOT system",
		10, -0.5, 9.5, 300, 0, 50},
	h1_foot_q{"FOOT Track Q [e]@Only SCI cuts, full data", 0xFF7C0A_c, cfg.all_q_binning},
	h1_track_q{"Track charge [Q]@Full FOOT system", HistConfig::fill, cfg.q_binning },
	h2_track_distance{"Track distances to vertex [mm]:Track Index [0,1,..N-1]@smaller x-axis means larger Q particle",
		track_count ? track_count : 10,
		-0.5,
		track_count ? track_count-0.5 : 9.5,
		600, 0, std::isfinite(cfg.distance_cut) ? (2*cfg.distance_cut): 30.0},
	h1_track_distance{"Track distances to common vertex [mm]", 0x99004c_c,
		200, 0, std::isfinite(cfg.distance_cut) ? (2*cfg.distance_cut): 5.0},
	h1_angle_ex{"((h1_angle_ex)) #rho [mrad]@#rho value", HistConfig::fill,
		cfg.rho_binning[0], cfg.rho_binning[1], cfg.rho_binning[2]
	},
	h2_vertex_z{"#rho [mrad]:Vertex z [mm]@Traced by the FOOT", 140, -80, 80,
		cfg.rho_binning[0], cfg.rho_binning[1], cfg.rho_binning[2]},
	h2_vertex_x{"#rho [mrad]:Vertex x [mm]@Traced by the FOOT", 140, -60, 60,
		cfg.rho_binning[0], cfg.rho_binning[1], cfg.rho_binning[2]},
	h2_nlayers_hit{"N layers hit in a track:Track Index [0,1,..N-1]@Full FOOT system, smaller x-axis means larger Q particle",
		10, -0.5, 9.5,   n_pairs, 0.5, n_pairs+0.5 },
	h2_rho_vs_theta{"#rho [mrad]:#theta(p,HI) [mrad]@Traced by the FOOT",
		cfg.theta_binning[0], cfg.theta_binning[1], cfg.theta_binning[2],
		cfg.rho_binning[0], cfg.rho_binning[1], cfg.rho_binning[2]} {
	if(std::isfinite(cfg.distance_cut) && track_count) h1_track_mult.AppendToTitle(" 0 means event not viable.");
}


Analysis::Analysis(const transport::Config& geometry,
	const HistConfig& cfg, const std::vector<Resonance>& channels, double nominal_beta_gamma)
	: h(cfg, std::isfinite(cfg.distance_cut) ? common_track_count(channels) : 0),
	  config(cfg), beta_gamma(nominal_beta_gamma) {
	if(std::any_of(channels.begin(), channels.end(), [](const auto& r) { return r.daughters.size() == 3; }))
		h.theta.emplace(cfg);

	if(cfg.do_scaling) {
		h.h1_angle_ex.AppendToAxisTitle(" (beta corrected)");
		for(auto* hist : {&h.h2_vertex_x, &h.h2_vertex_z, &h.h2_rho_vs_theta}) {
			auto* axis = (*hist)->GetYaxis();
			axis->SetTitle(Form("%s (beta corrected)", axis->GetTitle()));
		}
		if(h.theta) {
			auto* axis = h.theta->rho_vs_opening->GetYaxis();
			axis->SetTitle(Form("%s (beta corrected)", axis->GetTitle()));
		}
	}

	Eigen::Matrix<double, 8, 4> design;
	
	for(size_t i = 0; i < geometry.planes.size(); ++i) {
		const auto& p = geometry.planes[i];
		design.row(i) << p.a, p.a*p.z_mm, p.b, p.b*p.z_mm;
	}

	auto qr = design.colPivHouseholderQr();
	if(qr.rank() != 4)
		ERROR("Degenerate FOOT telescope geometry?");
	
	fit = qr.solve(Eigen::Matrix<double, 8, 8>::Identity());
}

/* Take the event and sort out its values into different histo's */
void Analysis::record(const ReactionEvent& e) {
	double rho = NAN;
	std::vector<Line3D> lines;
	std::vector<unsigned> charges;
	std::array<bool, 8> plane_fired{};
	for(size_t i = 0; i < e.measurements.size(); ++i) {
		const auto& hits = e.measurements[i];
		for(size_t j = 0; j < hits.size(); ++j)
			plane_fired[j] |= std::isfinite(hits[j]);
		Eigen::Map<const Eigen::Matrix<double, 8, 1>> measurement(hits.data());

		if(!measurement.allFinite()) continue;

		const Eigen::Vector4d t = fit * measurement;
		lines.emplace_back(t[0], t[1], t[2], t[3]);
		charges.push_back(e.channel->daughters[i].Z);
	}
	int pairs = 0;
	for(size_t i = 0; i < plane_fired.size(); i += 2) pairs += plane_fired[i] && plane_fired[i+1];
	h.h1_foots_fired.Fill(pairs);
	const auto track_count = lines.size();
	for(auto q: charges) {
		h.h1_foot_q.Fill(q);
	}
	if(e.channel && lines.size() == e.channel->daughters.size()) {
		const auto vertex = mnd::geom::FindVertex(lines);
		if(!vertex.is_null() && (!std::isfinite(config.distance_cut) ||
			std::all_of(lines.begin(), lines.end(), [&](const auto& line) {
				return line.DistanceTo(vertex) <= config.distance_cut;
			}))) {
			rho = phy::rho(phy::RhoExpressionType::p, e.channel->daughters, lines);
			if(config.do_scaling)
				rho *= e.mother_beta_gamma / beta_gamma;
			h.h1_angle_ex.Fill(rho);
			h.h2_vertex_z.Fill(vertex.z, rho);
			h.h2_vertex_x.Fill(vertex.x, rho);
			std::vector<double> proton_angles;
			for(size_t i = 0; i < lines.size(); ++i) {
				const double d = lines[i].DistanceTo(vertex);
				h.h1_track_distance.Fill(d);
				h.h2_track_distance.Fill(i, d);
				h.h2_nlayers_hit.Fill(i, n_pairs);
				h.h1_track_q.Fill(charges[i]);
				h.h2_q_vs_mult.Fill(track_count, charges[i]);
				if(i) {
					proton_angles.push_back(phy::MRAD_CVT * lines[i].AngleRelativeTo(lines.front()));
					h.h2_rho_vs_theta.Fill(proton_angles.back(), rho);
				}
			}
			if(h.theta && proton_angles.size() == 2) {
				h.theta->theta1_vs_theta2.Fill(proton_angles[1], proton_angles[0]);
				h.theta->rho_vs_opening.Fill(phy::MRAD_CVT * lines[1].AngleRelativeTo(lines[2]), rho);
			}
		}
	}
	h.h1_track_mult.Fill(std::isfinite(config.distance_cut) && !std::isfinite(rho) ? 0 : track_count);
}

void Analysis::write() {
	TCanvas cm("Multp", "Different tracks and extra info", 2150, 1400);
	cm.Divide(3,2);
	cm.cd(1); h.h2_q_vs_mult.Draw("COLZ"); gPad->SetLogz();
	cm.cd(2); h.h2_track_distance.Draw("COLZ"); gPad->SetLogz();
	// Truth-associated tracks have no reconstruction score.
	cm.cd(3); h.h2_score_vs_mult.Draw("COLZ"); gPad->SetLogz();
	cm.cd(4); h.h1_track_mult.Draw(); gPad->SetLogy();
	cm.cd(5); h.h2_nlayers_hit.Draw("COLZ"); gPad->SetLogz();
	cm.cd(6); h.h1_foots_fired.Draw(); gPad->SetLogy();

	const double rho_to_mev = phy::mp * phy::nuc::c * phy::nuc::c / 2 *
		beta_gamma * beta_gamma / (phy::MRAD_CVT * phy::MRAD_CVT);
	const auto to_energy = [rho_to_mev](double angle) { return rho_to_mev * angle * angle; };
	const auto to_angle = [rho_to_mev](double energy) { return sqrt(energy / rho_to_mev); };
	TCanvas ct("Physics", "Recognized tracks and angles", 2150, 1400);
	ct.Divide(2,2);
	ct.cd(1); h.h2_rho_vs_theta.Draw("COLZ");
	ct.cd(2); h.h2_vertex_x.Draw<TH2P::Y>(to_energy, to_angle, "Q_{#perp}  [MeV]", "COLZ");
	ct.cd(3); h.h1_angle_ex.Draw(to_energy, to_angle, "Q_{#perp}  [MeV]", "HIST");
	h.h1_angle_ex->SetStats(std::get<2>(config.rho_style));
	gPad->Update();
	if(auto* title = dynamic_cast<TPaveText*>(gPad->GetPrimitive("title")); title) {
		title->SetTextFont(std::get<0>(config.rho_style));
		title->SetTextSize(std::get<1>(config.rho_style));
		title->SetBorderSize(6);
		title->SetShadowColor((0x343434_c).GetColorCode());
	}
	ct.cd(4); h.h2_vertex_z.Draw<TH2P::Y>(to_energy, to_angle, "Q_{#perp}  [MeV]", "COLZ");

	TCanvas ct2("Physics2", "Recognized tracks and more", 2150, 1400);
	ct2.Divide(2,2);
	ct2.cd(1); h.h1_foot_q.Draw();
	ct2.cd(2); h.h2_track_distance.Draw("COLZ");
	ct2.cd(3); h.h1_track_q.Draw();
	ct2.cd(4); h.h1_track_distance.Draw();

	std::optional<TCanvas> ct3;
	if(h.theta) {
		ct3.emplace("Physics3", "Angles correction", 1900, 1600);
		ct3->Divide(2,1);
		ct3->cd(1); h.theta->theta1_vs_theta2.Draw("COLZ");
		ct3->cd(2); h.theta->rho_vs_opening.Draw("COLZ");
	}

	canvas::save_all<canvas::Exe>(config.save, mnd::to_views(config.info), { "simulated"sv } );
	// The histograms own the axis mappings. ROOT's global list only needs them while drawing.
	for(auto* hist : std::array<TH1*, 3>{&h.h2_vertex_x.h, &h.h1_angle_ex.h, &h.h2_vertex_z.h}) {
		for(auto* obj : *hist->GetListOfFunctions()) {
			if(auto* mapping = dynamic_cast<TF1*>(obj)) mapping->AddToGlobalList(false);
		}
	}
}
