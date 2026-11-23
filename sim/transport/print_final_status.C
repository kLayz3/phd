#include "../../includes/monad/accessories/formatting.hxx"
#include "../../includes/monad/accessories/meta.hxx"
#include "../../includes/nlohmann/json.hpp"

#include <TFile.h>
#include <TNamed.h>
#include <array>
#include <memory>
#include <string>

// root -l -b -q 'sim/transport/print_final_status.C("sim/gen/analyse_transport/9C_trial.root")'
void print_final_status(const char* filename = "") {
	/* Compare fitted residual widths with the foot_cov.json reference snapshot
	 * stored by analyse_transport. Trial sigmas are intrinsic noise inputs;
	 * residuals include both transport and track-prediction uncertainty. */
	if(!filename || !*filename)
		ERROR("Pass the final analyse_transport ROOT file");

	std::unique_ptr<TFile> file{TFile::Open(filename, "READ")};
	if(!file || file->IsZombie())
		ERROR("Cannot open {}", filename);
	auto* description = file->Get<TNamed>("analysis_summary");
	if(!description)
		ERROR("Missing analysis_summary in {}", filename);

	const auto summary = nlohmann::json::parse(description->GetTitle());
	const auto& sigmas = summary.at("sigma_um");
	const auto& residuals = summary.at("residuals");
	if(!sigmas.is_array() || sigmas.size() != 8 ||
		!residuals.is_array() || residuals.size() != 8)
		ERROR("Expected eight trial sigmas and eight residuals in {}", filename);

	std::array<const nlohmann::json*, 8> rows{};
	for(const auto& residual : residuals) {
		const u32 detector = residual.at("detector").get<u32>();
		const u32 pair = residual.at("pair").get<u32>();
		const auto axis = residual.at("axis").get<std::string>();
		if(detector >= rows.size() || pair >= 4 || rows[detector] ||
			(axis != "u" && axis != "v"))
			ERROR("Invalid or duplicate residual detector mapping in {}", filename);
		rows[detector] = &residual;
	}

	fmt::println("File: {}\nTransport input: {}\nSection: {}", filename,
		summary.at("input").get<std::string>(), summary.at("section").get<std::string>());
	fmt::println("Widths and differences in µm; delta = trial residual - measured residual.");
	fmt::println("Measured values are the reference snapshot saved with this trial.\n");
	fmt::println("{:>8} {:>4} {:>4} {:>12} {:>14} {:>14} {:>12} {:>10}",
		"Detector", "Pair", "Axis", "Input sigma", "Trial residual", "Measured", "Delta", "Delta [%]");

	for(size_t i = 0; i < rows.size(); ++i) {
		const auto& row = *rows[i];
		const double sigma = sigmas.at(i).get<double>();
		if(!mnd::is_non_negative(sigma))
			ERROR("Invalid trial sigma for detector d{}", i);

		auto width_um = [&](const char* key) {
			const auto& value = row.at(key);
			return value.is_number() ? value.get<double>() * 1000 : NAN;
		};
		const double trial = width_um("core_sigma_mm");
		const double measured = width_um("reference_core_sigma_mm");
		const bool fitted = mnd::is_positive(trial);
		const bool referenced = mnd::is_positive(measured);
		const std::string delta = fitted && referenced
			? fmt::format("{:+.2f}", trial - measured) : "n/a";
		const std::string percent = fitted && referenced
			? fmt::format("{:+.2f}", 100 * (trial - measured) / measured) : "n/a";

		fmt::println("{:>8} {:>4} {:>4} {:>12.2f} {:>14} {:>14} {:>12} {:>10}",
			fmt::format("d{}", i), row.at("pair").get<u32>() + 1,
			row.at("axis").get<std::string>(), sigma,
			fitted ? fmt::format("{:.2f}", trial) : "fit failed",
			referenced ? fmt::format("{:.2f}", measured) : "n/a", delta, percent);
	}
}
