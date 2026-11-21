#include "RtypesCore.h"

#include "monad/accessories/timing.hxx"
#include "monad/accessories/meta.hxx"

#include "util/CLI.h"
#include "util/TextUtils.h"
#include "util/FilesystemUtils.h"
#include "util/GaussFitMax.hxx"
#include "nlohmann/json.hpp"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TTree.h"
#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>

/* Okay, so this part uses no Geant4 anymore. All of the Geant4 stuff, related to straggling
 * must have been already precomputed by the transport program. Here we just "analyse it" ...
 * Basically trying to */
int main(int argc, char** argv) {
	std::string input, output, pdf, config = FOOT_TRACK_CONFIG, section = "c9";
	
	std::array<double, 3> bins_u{100, -5, 5}, bins_v{100, -5, 5};
	double side = 2.0;
	bool no_config = false;
	std::vector<double> sigma_um(8, 0);
	
	CLI::App app{"Four-pair leave-one-pair-out residuals from FOOT transport truth"};
	
	add_logged_option(app, "-f,--file", input, "Transport ROOT input (requires its .json summary)")
		->required()
		->check(CLI::ExistingFile);
	add_logged_option(app, "-o,--output", output, "Residual ROOT output; default INPUT_residuals.root");
	add_logged_option(app, "--pdf", pdf, "Residual canvas PDF; default derived from output");
	add_logged_option(app, "-c,--config", config, "foot_track_analysis configuration");
	add_logged_option(app, "--section", section, "Configuration section (actual carbon section is c9)");
	add_logged_flag(app, "--no-config", no_config, "Use the compiled foot_track_analysis defaults");
	auto* u_option = add_logged_option(app, "--binning-dx", bins_u, "Residual u: bins,min,max [mm]")
		->delimiter(',');
	auto* v_option = add_logged_option(app, "--binning-dy", bins_v, "Residual v: bins,min,max [mm]")
		->delimiter(',');
	auto* side_option = add_logged_option(app, "-s,--side-ratio", side, "GaussFitMax core-fit side ratio")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--sigma-um", sigma_um, "Trial detector noise RMS, d0..d7 [micrometres]; default zero")
		->delimiter(',')
		->expected(8);
	
	CLI11_PARSE(app, argc, argv);
	srand(time(NULL));

	std::vector<TimePoint> tv;
	tv.emplace_back(TimePoint("start"));
	
	u64 seed = (u64)rand() | (u64)rand();

	// Match quick/run: compiled defaults, GLOBAL, SECTION; explicit analysis arguments override these.
	const auto cli_u = bins_u, cli_v = bins_v;
	const auto cli_side = side;
	if(!no_config) {
		const auto text = mnd::fs::parse_file_to_string(config);
		std::string body{mnd::sdsl::extract_text_body("GLOBAL", "", text).value_or("")};
		if(!section.empty()) {
			const auto selected = mnd::sdsl::extract_text_body("SECTION", section, text);
			if(selected.is_none()) ERROR("Unknown config section: {}", section);
			body += "\n" + std::string(selected.unwrap());
		}
		CLI::App settings;
		settings.allow_extras(); // Input files and charge gates belong to the real-data analysis.
		add_logged_option(settings, "--binning-dx", bins_u, "")
			->delimiter(',');
		add_logged_option(settings, "--binning-dy", bins_v, "")
			->delimiter(',');
		add_logged_option(settings, "-s,--side-ratio", side, "");
		auto args = mnd::parse_argv(body, "foot_track_analysis settings");
		settings.parse(args.argc(), args.data());
	}
	if(u_option->count()) bins_u = cli_u;
	if(v_option->count()) bins_v = cli_v;
	if(side_option->count()) side = cli_side;
	for(const auto& bins: {bins_u, bins_v})
		if(!mnd::isfinite(bins[0], bins[1], bins[2]) || bins[0] < 1 || bins[0] > 1000000 || std::floor(bins[0]) != bins[0] ||
			bins[1] >= bins[2])
			ERROR("Invalid residual binning");
	if(!mnd::is_positive(side)) ERROR("Invalid core-fit side ratio");
	for(double sigma: sigma_um)
		if(!mnd::is_non_negative(sigma)) ERROR("Noise RMS must be finite and nonnegative");
	if(output.empty()) output = (std::filesystem::path(input).parent_path() /
		(std::filesystem::path(input).stem().string()+"_residuals.root")).string();
	if(pdf.empty()) {
		auto p = std::filesystem::path(output);
		p.replace_extension(".pdf");
		pdf = p.string();
	}
	if(output == pdf || std::filesystem::exists(output) || std::filesystem::exists(pdf))
		ERROR("Output already exists or ROOT/PDF paths coincide");
	nlohmann::json metadata;
	std::ifstream info(input+".json");
	if(!info) ERROR("Cannot read transport metadata: {}.json", input);
	info >> metadata;
	using Matrix = Eigen::Matrix<double, 8, 8>;
	Eigen::Matrix<double, 8, 4> A;
	Eigen::Matrix<double, 8, 1> offsets;
	std::array<u32, 8> pairs;
	std::array<bool, 8> is_u, seen{};
	std::array<std::array<bool, 2>, 4> axes{};
	if(metadata.at("planes").size() != 8) ERROR("Expected eight plane descriptions");
	for(const auto& plane: metadata.at("planes")) {
		const u32 i = plane.at("index"), pair = plane.at("pair");
		const std::string axis = plane.at("axis");
		if(i >= 8 || pair >= 4 || seen[i] || (axis != "u" && axis != "v"))
			ERROR("Invalid plane mapping in metadata");
		seen[i] = true;
		pairs[i] = pair;
		is_u[i] = axis == "u";
		if(axes[pair][is_u[i]]) ERROR("Duplicate pair readout axis");
		axes[pair][is_u[i]] = true;
		const double a = plane.at("a"), b = plane.at("b"), z = plane.at("z_mm");
		A.row(i) << a, a*z, b, b*z;
		offsets(i) = plane.at("offset_mm");
	}
	if(!A.allFinite() || !offsets.allFinite()) ERROR("Nonfinite geometry");
	Matrix R = Matrix::Identity();
	for(u32 pair=0; pair<4; ++pair) {
		Eigen::Matrix<double, 6, 4> reduced;
		std::array<u32, 6> kept;
		u32 row = 0;
		for(u32 i=0; i<8; ++i) if(pairs[i] != pair) {
			kept[row] = i;
			reduced.row(row++) = A.row(i);
		}
		auto qr = reduced.colPivHouseholderQr();
		if(qr.rank() != 4) ERROR("Rank-deficient leave-one-pair-out fit");
		const Eigen::Matrix<double, 4, 6> L = qr.solve(Eigen::Matrix<double, 6, 6>::Identity());
		for(u32 i=0; i<8; ++i) if(pairs[i] == pair)
			for(u32 j=0; j<6; ++j) R(i, kept[j]) = -(A.row(i)*L.col(j))(0);
	}
	const Matrix T = R.array().square().matrix();
	Eigen::JacobiSVD<Matrix> svd(T);
	const auto singular = svd.singularValues();
	const u32 rank = (singular.array() > singular(0)*1e-10).count();
	WARN("Effective residual binning u: {}; v: {}; side ratio: {:.6g}\n"
		"Variance-response rank (relative tolerance 1e-10): {}/8\n", bins_u, bins_v, side, rank);
	if(rank < 8) WARN("Eight independent detector variances are not identified at this tolerance.\n");

	std::unique_ptr<TFile> file{TFile::Open(input.c_str(), "READ")};
	if(!file || file->IsZombie()) ERROR("Cannot open transport ROOT input");
	auto* tree = file->Get<TTree>("crossings");
	if(!tree) ERROR("Missing crossings tree");
	std::array<int, 8> present{};
	std::array<double, 8> measurement{};
	std::array<std::unique_ptr<TH1D>, 8> hist;
	tree->SetBranchStatus("*", 0);
	for(u32 i=0; i<8; ++i) {
		const auto prefix = "d"+std::to_string(i)+"_";
		for(const auto& name: {prefix+"present", prefix+"measurement_mm"}) tree->SetBranchStatus(name.c_str(), 1);
		if(tree->SetBranchAddress((prefix+"present").c_str(), &present[i]) < 0 ||
			tree->SetBranchAddress((prefix+"measurement_mm").c_str(), &measurement[i]) < 0)
			ERROR("Invalid crossing branch");
		const auto& bins = is_u[i] ? bins_u : bins_v;
		const auto name = "residual_"+prefix+(is_u[i] ? "u" : "v");
		const auto title = "Pair "+std::to_string(pairs[i])+" "+(is_u[i] ? "u" : "v")+
			" (d"+std::to_string(i)+");Measurement - prediction [mm];Count";
		hist[i] = std::make_unique<TH1D>(name.c_str(), title.c_str(), static_cast<int>(bins[0]), bins[1], bins[2]);
		hist[i]->SetDirectory(nullptr);
		hist[i]->SetStatOverflows(TH1::kConsider);
		hist[i]->SetFillColor(kRed-9);
		hist[i]->SetFillStyle(1001);
	}
	std::mt19937_64 random(seed);
	std::normal_distribution<double> normal;
	Long64_t accepted = 0;
	Eigen::Matrix<double, 8, 1> mean = Eigen::Matrix<double, 8, 1>::Zero();
	Matrix sum_cov = Matrix::Zero();

	const Long64_t nentries = tree->GetEntries();
	for(Long64_t event=0; event < nentries; ++event) {
		if(tree->GetEntry(event) <= 0) ERROR("Cannot read transport entry");
		bool valid = true;
		Eigen::Matrix<double, 8, 1> measured;
		for(u32 i=0; i<8; ++i) {
			valid &= present[i] && mnd::isfinite(measurement[i]);
			measured(i) = measurement[i]+offsets(i);
		}
		if(!valid) continue;
		for(u32 i=0; i<8; ++i) measured(i) += normal(random)*sigma_um[i]*0.001;
		const Eigen::Matrix<double, 8, 1> residual = R*measured;
		++accepted;
		const Eigen::Matrix<double, 8, 1> delta = residual-mean;
		mean += delta/static_cast<double>(accepted);
		sum_cov += delta*(residual-mean).transpose();
		for(u32 i=0; i<8; ++i) hist[i]->Fill(residual(i));
	}
	if(!accepted) ERROR("No finite four-pair crossings available");
	gROOT->SetBatch(true);
	
	TCanvas canvas("TransportResiduals", "Four-pair leave-one-pair-out residuals", 1800, 950);
	canvas.Divide(4, 2);
	
	nlohmann::json summary {
		{"input",                  input},
		{"config",                 no_config ? "" : config},
		{"section",                section},
		{"binning_dx",             bins_u},
		{"binning_dy",             bins_v},
		{"side_ratio",             side},
		{"fit_iterations",         4},
		{"sigma_um",               sigma_um},
		{"seed",                   seed},
		{"four_pair_tracks",       accepted},
		{"variance_response_rank", rank},
		{"selection", "All finite four-pair transport events; no real-data charge gates, CKF pruning or final score cut"}
	};
	TFile destination(output.c_str(), "RECREATE");
	if(destination.IsZombie()) ERROR("Cannot create residual output");
	for(u32 i=0; i<8; ++i) {
		canvas.cd(pairs[i]+1+(is_u[i] ? 0 : 4));
		hist[i]->Draw("HIST");
		const auto fit = GaussFitMax(hist[i].get(), side, 4);
		const bool fitted = VerifyGaussFitMax(fit);
		if(fitted) {
			TF1 gaussian((std::string(hist[i]->GetName())+"_core").c_str(), "gaus",
				fit.first[1]-side*fit.first[2], fit.first[1]+side*fit.first[2]);
			gaussian.SetParameters(fit.first.data());
			gaussian.DrawCopy("SAME");
		}
		double reference = NAN;
		for(const auto& item: metadata.at("reference_residuals"))
			if(item.at("pair").get<u32>() == pairs[i] && item.value("available", true))
				reference = item.at(is_u[i] ? "u_mm" : "v_mm");
		TLatex text;
		text.SetNDC();
		text.SetTextSize(.045);
		const std::string label = fmt::format("Core #sigma: {:.2f} #mum; measured: {:.2f} #mum",
			fitted ? fit.first[2]*1000 : NAN, reference*1000);
		text.DrawLatex(.12, .83, label.c_str());
		summary["residuals"].push_back({
			{"detector",                i},
			{"pair",                    pairs[i]},
			{"axis",                    is_u[i] ? "u" : "v"},
			{"mean_mm",                 mean(i)},
			{"variance_mm2",            accepted>1 ? sum_cov(i, i)/(accepted-1) : NAN},
			{"core_sigma_mm",           fitted ? fit.first[2] : NAN},
			{"reference_core_sigma_mm", reference}
		});
		WARN("{}: core={:.6g} um; reference={:.6g} um\n", hist[i]->GetName(),
			fitted ? fit.first[2]*1000 : NAN, reference*1000);
		hist[i]->Write();
	}
	auto write_matrix = [&](const char* name, const Matrix& value) {
		TMatrixD matrix(8, 8);
		for(u32 i=0; i<8; ++i) for(u32 j=0; j<8; ++j) matrix(i, j)=value(i, j);
		matrix.Write(name);
	};
	write_matrix("residual_response", R);
	write_matrix("variance_response", T);
	if(accepted>1) write_matrix("residual_covariance_mm2", sum_cov/(accepted-1));
	TNamed description("analysis_summary", summary.dump(4).c_str());
	description.Write();
	canvas.Write();
	canvas.Print(pdf.c_str());
	WARN("Used {}/{} events; wrote {} and {}\n", accepted, tree->GetEntries(), output, pdf);

	tv.emplace_back(TimePoint("end"));
	PrintElapsed<kSECOND>(std::move(tv));
}
