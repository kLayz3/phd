#include "util/RunsheetParser.h"
#include "util/JSONParser.h"
#include "util/FromChars.h"
#include <string>

nlohmann::json mnd::fs::runsheet_obj {};
constexpr static const char* runsheet_file_path_ = "params/runsheet.json";

/* Returns None for missing key, null or JSON conversion exception.
 * Otherwise return (by value) the json object converted to `T`. */
template<typename T, typename Key>
static mnd::Option<T> j_value(const nlohmann::json& j, const Key& key) {
	const auto it = j.find(key);
	if(it == j.end() || it->is_null())
		return mnd::None;
	try {
		return mnd::Some{ it->template get<T>() };
	} catch(const nlohmann::json::exception& ) {
		return mnd::None;
	}
}

/* Will throw on bad parse, or if file does not exist. */
void mnd::fs::load_runsheet(const std::filesystem::path& p) {
	runsheet_obj = ParseJSON(p);
}
void mnd::fs::load_runsheet() {
	runsheet_obj = ParseJSON(runsheet_file_path_);
}

static const auto is_eq = [](double a, double b) -> bool {
	return a == b || (std::isnan(a) && std::isnan(b));
};

bool mnd::RunsheetState::Brho::operator==(const Brho& rhs) const noexcept {
	return is_eq(s1_s2, rhs.s1_s2) &&
	       is_eq(s2_s3, rhs.s2_s3) &&
	       is_eq(s3_s4, rhs.s3_s4);
}
bool mnd::RunsheetState::operator==(const RunsheetState& rhs) const noexcept {
	return brho == rhs.brho &&
	       is_eq(e0, rhs.e0) &&
	       is_eq(primary.Z, rhs.primary.Z) &&
	       is_eq(primary.A, rhs.primary.A) &&
	       is_eq(secondary.Z, rhs.secondary.Z) &&
	       is_eq(secondary.A, rhs.secondary.A) &&
	       secondary.N_electrons == rhs.secondary.N_electrons;
}
bool mnd::RunsheetState::operator!=(const RunsheetState& rhs) const noexcept {
	return !(*this == rhs);
}

mnd::RunsheetState mnd::RunsheetState::from(const nlohmann::json& j) {
	return RunsheetState {
		.brho = {
			.s1_s2 = j_value<double>(j, fs::brho::s1_s2).value_or(NAN),
			.s2_s3 = j_value<double>(j, fs::brho::s2_s3).value_or(NAN),
			.s3_s4 = j_value<double>(j, fs::brho::s3_s4).value_or(NAN)
		},
		.e0 = j_value<double>(j, fs::e_primary).value_or(NAN),
		.primary  = phy::Nucleus::get_ion(
			j_value<std::string>(j, fs::i_primary).value_or(""), false, fs::i_primary
		),
		.secondary = phy::Nucleus::get_ion(
			j_value<std::string>(j, fs::i_secondary).value_or(""), false, fs::i_secondary
		)
	};
}

std::ostream& mnd::operator<<(std::ostream& os, const RunsheetState& r) {
	return os << KBH_GRN "({" KRNM " brho: {"
		<< fs::brho::s1_s2 << ": " << BOLD << r.brho.s1_s2 << KNRM " Tm, "
		<< fs::brho::s2_s3 << ": " << BOLD << r.brho.s2_s3 << KNRM " Tm, "
		<< fs::brho::s3_s4 << ": " << BOLD << r.brho.s3_s4 << KNRM " Tm}"
		<< ", {primary"
		<< ": Z: " BOLD << r.primary.Z << KNRM
		<< ", A: " BOLD << r.primary.A << KNRM
		<< ", Q: " BOLD << r.primary.charge_state().map([](auto x) { return std::to_string(x); }).value_or("none")
		<< KNRM "[e], EKin: " << KBH_RED << r.e0 << KNRM " MeV/u}"
		<< ", {secondary"
		<< ": Z: " BOLD << r.secondary.Z << KNRM
		<< ", A: " BOLD << r.secondary.A << KNRM
		<< ", Q: " BOLD << r.secondary.charge_state().map([](auto x) { return std::to_string(x); }).value_or("none")
		<< KNRM "[e]} " KBH_GRN "})" KNRM;
}

template<>
mnd::OptRunsheetStatePair
mnd::QueryRunsheet<
	false,
	mnd::OptRunsheetStatePair
> (std::string_view str) {
	if(fs::runsheet_obj.empty())
		MND_THROW("Runsheet object is empty (null). Did you load the runsheet?");

	const auto [
		basename,
		start_num,
		end_num
	] = mnd::fs::file_info(str);
	
	Option<RunsheetState> start_info, end_info;

	/* Loop over the JSON object, select the rows named as the filename
	 * base and look from there. */
	for(const auto& [_, row] : fs::runsheet_obj.items()) {
		if(row.at(fs::file_name_key) != basename)
			continue;

		if(row.at(fs::start_num_key) == start_num)
			start_info = Some{ RunsheetState::from(row) };
		
		if(row.at(fs::end_num_key) == end_num)
			end_info = Some{ RunsheetState::from(row) };
	}

	return {
		std::move(start_info),
		std::move(end_info)
	};
}

template<>
mnd::RunsheetState mnd::QueryRunsheet<
	true,
	mnd::RunsheetState
> (std::string_view str) {
	/* Only used for debugging. */
	const auto [
		basename,
		start_num,
		end_num
	] = mnd::fs::file_info(str);

	auto [start_info, end_info] = QueryRunsheet<false, OptRunsheetStatePair>(str);

	if(start_info.is_none())
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename '%s' "
			"and start run-number '%u' not found in a row in the runsheet.\n",
			(int)str.size(), str.data(), basename.c_str(), start_num);

	if(end_info.is_none())
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename '%s' "
			"and end run-number '%u' not found in a row in the runsheet.\n",
			(int)str.size(), str.data(), basename.c_str(), end_num);

	if(start_info.unwrap() != end_info.unwrap()) {
		std::cerr << "mnd::QueryRunsheet: mismatch!\n"
			<< start_info.unwrap() << " ... and:\n" << end_info.unwrap() << std::endl;
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename '%s': "
			"start '%u' and end '%u' run-numbers mismatched runsheet info?\n",
			(int)str.size(), str.data(), basename.c_str(), start_num, end_num);
	}
	
	return start_info.unwrap();
}
