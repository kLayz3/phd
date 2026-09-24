#pragma once

#include <ostream>
#include <istream>
#include <string_view>
#include <filesystem>

#include "util/MacroHelpers.h"
#include "util/Option.hxx"
#include "util/MPhysics.h"

#include "nlohmann/json.hpp"

namespace mnd {

constexpr static u32 DEFAULT_Z_PRIMARY   =  6;
constexpr static u32 DEFAULT_A_PRIMARY   = 12;
constexpr static u32 DEFAULT_Z_SECONDARY =  6;
constexpr static u32 DEFAULT_A_SECONDARY =  9;

namespace fs {

inline constexpr const char* runsheet_file_path = "params/runsheet.json";
inline constexpr const char* file_name_key = "Name";
inline constexpr const char* start_num_key = "Start file number";
inline constexpr const char* end_num_key   = "Stop file number";
inline constexpr const char* e_primary     = "E_in [MeV/u]";
inline constexpr const char* i_primary     = "Ion";
inline constexpr const char* i_secondary   = "Fragment";

namespace brho {

inline constexpr const char* s1_s2 = "S1-S2";
inline constexpr const char* s2_s3 = "S2-S3";
inline constexpr const char* s3_s4 = "S3-S4";

} // namespace brho

/* Runsheet object won't ever be mutated. Safe to keep it non-thread_local.
 * It will only get loaded once inside the main function. */
extern ::nlohmann::json runsheet_obj;

/* This should be called exclusively by the main thread, only once. */
void load_runsheet(const std::string_view = runsheet_file_path);

} // namespace fs

/* Different interesting things we can query from a single runsheet row. */
struct RunsheetState {
	struct Brho {
		double s1_s2, s2_s3, s3_s4;
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
std::ostream& operator<<(std::ostream& , const RunsheetState& );

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
