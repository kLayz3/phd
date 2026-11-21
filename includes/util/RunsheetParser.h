#pragma once

#include <istream>
#include <string_view>
#include <filesystem>

#include "util/MacroHelpers.h"
#include "util/Option.hxx"
#include "util/MPhysics.h"

#include "nlohmann/json.hpp"

namespace mnd {

inline constexpr u32 DEFAULT_Z_PRIMARY   =  6;
inline constexpr u32 DEFAULT_A_PRIMARY   = 12;
inline constexpr u32 DEFAULT_Z_SECONDARY =  6;
inline constexpr u32 DEFAULT_A_SECONDARY =  9;

namespace fs {

inline constexpr const char* runsheet_file_path = "params/runsheet.json";
inline constexpr const char* file_name_key = "Name";
inline constexpr const char* start_num_key = "Start file number";
inline constexpr const char* end_num_key   = "Stop file number";
inline constexpr const char* e_primary     = "E_in [MeV/u]";
inline constexpr const char* i_primary     = "Ion";
inline constexpr const char* i_secondary   = "Fragment";

namespace brho {

inline constexpr const char* label[] = {
	"TA-S1",
	"S1-S2",
	"S2-S3",
	"S3-S4"
};
inline constexpr size_t label_index(std::string_view name) {
	for(size_t i = 0; i < std::size(label); ++i) {
		if(name == label[i])
			return i;
	}
	throw "size_t mnd::fs::brho::label_index(std::string_view ): Unknown label provided.";
}

} // namespace brho

/* Runsheet object won't ever be mutated. Safe to keep it non-thread_local.
 * It will only get loaded once inside the main function. */
extern ::nlohmann::json runsheet_obj;

/* This should be called exclusively by the main thread, only once. */
void load_runsheet(const std::string_view = runsheet_file_path);

} // namespace fs

/* Different interesting things we can query from a single runsheet row. */
inline constexpr const char* focal_point[] = {
	"S1", "S2", "S3", "S4"
};
inline constexpr size_t N_FOCAL_PTS = std::size(focal_point);

static_assert(std::size(fs::brho::label) == N_FOCAL_PTS,
	"Compiling in $1 focal points, but there's $2 labels in the array. Not matching."
);
inline constexpr size_t TA_S1 = fs::brho::label_index("TA-S1");
inline constexpr size_t S1_S2 = fs::brho::label_index("S1-S2");
inline constexpr size_t S2_S3 = fs::brho::label_index("S2-S3");
inline constexpr size_t S3_S4 = fs::brho::label_index("S3-S4");

struct RunsheetState {
	struct Brho {
		double value[N_FOCAL_PTS];
		
		/* Can throw on a bad string query, only proper mnd::fs::brho::label queries allowed. */
		inline constexpr double& operator[](std::string_view name) {
			return value[ fs::brho::label_index(name) ];
		};
		inline constexpr double const& operator[](std::string_view name) const {
			return value[ fs::brho::label_index(name) ];
		};
		inline constexpr double& operator[](size_t I) noexcept { return value[I]; };
		inline constexpr double const& operator[](size_t I) const noexcept { return value[I]; };
		bool operator==(const Brho& rhs) const noexcept;
	} brho;

	/* Primary beam energy in [MeV/u]. Read from the file. */
	double e0;
	::phy::Nucleus primary;
	
	/* Secondary beam main nucleus. */
	::phy::Nucleus secondary;

	bool operator==(const RunsheetState& rhs) const noexcept;
	bool operator!=(const RunsheetState& rhs) const noexcept;

	static RunsheetState from(const nlohmann::json &);
};

using OptRunsheetStatePair = std::pair<
	Option<RunsheetState>,
	Option<RunsheetState>
>;

/* Query the JSON with complete file name, which can be multiple files concatenated. */
template<bool DoCheck = true, typename ResultType = RunsheetState>
ResultType QueryRunsheet(std::string_view );

/* Return runsheet status at the enclosed run numbers.
 * Will check that these statuses match.
 * In case they don't match, an exception is thrown. */
template<>
[[nodiscard]] RunsheetState QueryRunsheet<true, RunsheetState>(std::string_view );

/* Return pair of possible runsheet statuses at the enclosed run numbers. */
template<>
[[nodiscard]] OptRunsheetStatePair QueryRunsheet<false, OptRunsheetStatePair>(std::string_view );

} // namespace mnd

template<>
struct fmt::formatter<mnd::RunsheetState::Brho> {
	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

	template<typename Context>
	auto format(const mnd::RunsheetState::Brho& brho, Context& ctx) const {
		auto out = fmt::format_to(ctx.out(), "brho: {{");
		for(size_t i = 0; i < mnd::N_FOCAL_PTS; ++i) {
			out = fmt::format_to(out, "{}{}: {:.6g} Tm", i ? ", " : "", mnd::fs::brho::label[i],
				styled(brho.value[i], fmt::emphasis::bold));
		}
		return fmt::format_to(out, "}}");
	}
};

template<>
struct fmt::formatter<mnd::RunsheetState> {
	constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

	template<typename Context>
	auto format(const mnd::RunsheetState& r, Context& ctx) const {
		constexpr auto bold = fmt::emphasis::bold;
		constexpr auto green = fg(fmt::terminal_color::bright_green) | bold;
		return fmt::format_to(ctx.out(), "{}{}, {{primary: Z: {}, A: {}, Q: {}[e], EKin: {:.6g} MeV/u}}"
			", {{secondary: Z: {}, A: {}, Q: {}[e]}} {}",
			styled("({", green), r.brho,
			styled(r.primary.Z, bold), styled(r.primary.A, bold),
			styled(r.primary.charge_state().map([](auto x) { return std::to_string(x); }).value_or("none"), bold),
			styled(r.e0, mnd::style::bold_red),
			styled(r.secondary.Z, bold), styled(r.secondary.A, bold),
			styled(r.secondary.charge_state().map([](auto x) { return std::to_string(x); }).value_or("none"), bold),
			styled("})", green));
	}
};
