#include "monad/accessories/formatting.hxx"
#include "monad/accessories/meta.hxx"

#include "TransportConfig.hh"
#include "util/JSONParser.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace transport {

Config load_config(const std::string& setup_path, const std::string& cov_path, const std::string& ion) {
	const auto setup = ParseJSON(setup_path);
	const auto cov = ParseJSON(cov_path);
	Config result;
	result.provenance = {
		{"sources", {
			{"setup", std::filesystem::absolute(setup_path).lexically_normal().string()},
			{"covariance", std::filesystem::absolute(cov_path).lexically_normal().string()}
			}
		},
		{"ion", ion},
		{"units", {
			{"position", "mm"}, {"slopes", "dimensionless"}, {"energy", "MeV"}
			}
		},
		{"measurement", "primary crossing at silicon midplane; bare detector-axis coordinate"}
	};
	u32 silicon = 0;
	
	for(const auto& item: setup.at("DSSSD").at("composition")) {
		auto name = item.at("material").get<std::string>();
		if(name == "Si") name = "G4_Si";
		const double thickness = item.at("thickness").get<double>() * 1000;
		if(!mnd::is_positive(thickness))
			ERROR("Material thickness must be positive and finite (input unit: metres)");
		silicon += name == "G4_Si";
		result.materials.push_back({name, thickness});
		result.provenance["composition_per_detector"].push_back({
			{"material", name}, {"thickness_mm", thickness}
		});
	}
	if(silicon != 1) ERROR("Composition must contain exactly one Si/G4_Si measurement slab");
	
	const auto& box = setup.at("box");
	const auto positions = box.at("det_pos").get<std::vector<double>>();
	if(positions.size() != N)
		ERROR("Expected {} detector positions, received: {}", N, positions.size());
	
	std::array<bool, N> seen{};
	std::array<std::array<bool, 2>, N/2> axes{};

	for(const auto& item: setup.items()) {
		if(item.key().compare(0, 4, "FOOT") != 0) continue;
		const auto& detector = item.value();
		const u32 index = detector.at("N").get<u32>();
		if(index >= N || seen[index]) ERROR("Invalid or duplicate FOOT index");
		seen[index] = true;
		
		const auto orientation = detector.at("orientation").get<std::string>();
		const auto axis = orientation.empty() ? ' ' : orientation.back();
		
		if(axis != 'x' && axis != 'y') throw std::runtime_error("Unknown detector orientation");

		const bool u = axis == 'x';
		if(axes[index/2][u])
			ERROR("Pair has duplicate readout axes");
		
		axes[index/2][u] = true;
		
		const double angle = detector.at("delta_a").get<double>() * M_PI/180;
		const double z = positions[index] + detector.value("dz", 0.0) - box.at("target").at("dz").get<double>();
		const double offset = detector.at("delta_p").get<double>();

		if(!mnd::isfinite(z, angle, offset))
			ERROR("Detector geometry (z,angle,offset) must be finite");

		result.planes[index] = {item.key(), index/2, u, z,
			u ? std::cos(angle) : -std::sin(angle),
			u ? std::sin(angle) : std::cos(angle), offset};
	}
	if( !std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }) )
		ERROR("Missing FOOT geometry");
	
	for(size_t i=1; i<N; ++i)
		if(result.planes[i].z_mm <= result.planes[i-1].z_mm)
			ERROR("Detector indices must be in increasing z order");

	// Residual widths are comparison data, never the input detector smearing.
	nlohmann::json reference = nlohmann::json::array();

	for(u32 pair=0; pair<N/2; ++pair) {
		const auto& entries = cov.at("values").at("pair " + std::to_string(pair+1)).at("position");
		bool found = false;
		for(const auto& entry: entries) {
			if(entry.at("ion") != ion)
				continue;
			if(found)
				ERROR("Duplicate residual species in a pair");
			found = true;
			
			const auto& r = entry.at("residue");
			const double u = r.at("u").get<double>() * 0.001;
			const double v = r.at("v").get<double>() * 0.001;
			if(!mnd::is_positive(u, v))
				ERROR("Residual widths must be positive and finite (input unit: υm)");
			reference.push_back({
				{"pair", pair}, {"u_mm", u}, {"v_mm", v}
			});
		}

		if(!found)
			reference.push_back({{"pair", pair}, {"available", false}});
	}

	result.provenance["reference_residuals"] = reference;
	for(size_t i=0; i<N; ++i) {
		const auto& p = result.planes[i];
		result.provenance["planes"].push_back({
			{"index", i}, {"name", p.name}, {"pair", p.pair},
			{"axis", p.u ? "u" : "v"}, {"z_mm", p.z_mm}, {"a", p.a}, {"b", p.b}, {"offset_mm", p.offset_mm}
		});
	}
	return result;
}

} // namespace transport
