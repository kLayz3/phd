#pragma once
#include "Decay.hh"
#include "TransportConfig.hh"
#include "util/PrettyHisto.h"
#include "util/MacroHelpers.h"
#include <optional>
#include <tuple>
#include "Eigen/Core"

#include <string>

struct HistConfig {
	inline static constexpr auto fill = mnd::col::literals::operator""_c(0x91E965)
		+ mnd::col::literals::operator""_o(0.77L);
	A3 theta_binning{200, 0, 200}, rho_binning{300, 0, 300};
	A3 q_binning{80, 0, 7.5}, all_q_binning{700, 0, 7.5};
	std::tuple<u32, double, bool> rho_style{132, 0.045, false};
	double distance_cut = NAN;
	bool do_scaling = false;
	std::vector<canvas::Extension> save;
	std::vector<std::string> info;

	void add_options(CLI::App& app);
};

struct Histos {
	TH1P h1_track_mult;
	TH1P h1_foots_fired;
	TH2P h2_q_vs_mult;
	TH2P h2_score_vs_mult;
	TH1P h1_foot_q;
	TH1P h1_track_q;
	TH2P h2_track_distance;
	TH1P h1_track_distance;
	TH1P h1_angle_ex;
	TH2P h2_vertex_z;
	TH2P h2_vertex_x;
	TH2P h2_nlayers_hit;
	TH2P h2_rho_vs_theta;
	struct ThetaHistograms {
		TH2P theta1_vs_theta2, rho_vs_opening;
		explicit ThetaHistograms(const HistConfig& cfg);
	};
	std::optional<ThetaHistograms> theta;

	Histos(const HistConfig& cfg, size_t track_count);
};

struct Analysis {
	Eigen::Matrix<double, 4, 8> fit;
	
	Histos h;
	const HistConfig& config;
	double beta_gamma;
	
	Analysis(const transport::Config& geometry,
		const HistConfig& cfg, const std::vector<Resonance>& channels, double nominal_beta_gamma);

	void record(const ReactionEvent& event);
	void write();
};
