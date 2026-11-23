#pragma once

#include "monad/accessories/primitives.hxx"
#include "nlohmann/json.hpp"
#include <array>
#include <string>
#include <vector>

namespace transport {

inline constexpr size_t N = 8;

struct MaterialLayer {
	std::string name;
	double thickness_mm;
};
struct Plane {
	std::string name;
	u32 pair;
	bool u;
	double z_mm, a, b, offset_mm;
};
struct Config {
	std::array<Plane, N> planes;
	std::vector<MaterialLayer> materials;
	nlohmann::json provenance;
};

Config load_config(const std::string& setup, const std::string& ion);

} // namespace transport
