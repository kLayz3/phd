#pragma once

#include <ostream>
#include <string_view>

#include "util/MacroHelpers.h"
#include "util/Option.hxx"
#include "MPhysics.hxx"

#include "nlohmann/json.hpp"

namespace mnd {

constexpr static u32 DEFAULT_Z_PRIMARY   =  6;
constexpr static u32 DEFAULT_A_PRIMARY   = 12;
constexpr static u32 DEFAULT_Z_SECONDARY =  6;
constexpr static u32 DEFAULT_A_SECONDARY =  9;

namespace fs {

inline constexpr const char* file_name_key = "Name ";
inline constexpr const char* start_num_key = "Start file number";
inline constexpr const char* stop_num_key  = "Stop file number";
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
void load_runsheet(const std::string& );

} // namespace fs

/* Different interesting things we can query from a single runsheet row. */
struct RunsheetState {
	struct Brho {
		double s1_s2, s2_s3, s3_s4;
		inline bool operator==(const Brho& rhs) const noexcept {
			return s1_s2 == rhs.s1_s2 &&
			       s2_s3 == rhs.s2_s3 &&
			       s3_s4 == rhs.s3_s4;
		};
	} brho;

	/* Primary beam energy in [MeV/u]. Read from the file. */
	double e0;
	::phy::Nucleus primary;
	::phy::Nucleus secondary;

	/* Secondary beam */
	inline bool operator==(const RunsheetState& rhs) const noexcept {
		return brho == rhs.brho &&
		       e0   == rhs.e0 &&
			   primary.Z == rhs.primary.Z &&
			   primary.A == rhs.primary.A &&
			   secondary.Z == rhs.secondary.Z &&
			   secondary.A == rhs.secondary.A &&
			   secondary.N_electrons == rhs.secondary.N_electrons;
	}
	inline bool operator!=(const RunsheetState& rhs) const noexcept {
		return !(*this == rhs);
	}

	static RunsheetState from(const nlohmann::json &);
	static ::phy::Nucleus get_ion(const std::string& , const char* );
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
 * If they don't, an exception is thrown. */
template<> RunsheetState QueryRunsheet<true, RunsheetState>(std::string_view );

/* Return pair of possible runsheet statuses at the enclosed run numbers. */
template<> OptRunsheetStatePair QueryRunsheet<false, OptRunsheetStatePair>(std::string_view );

} // namespace mnd
