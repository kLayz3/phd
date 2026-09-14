#include "RunsheetParser.h"
#include "JSONParser.h"
#include "FromChars.h"
#include "MPhysics.hxx"
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

/* Parse an nucleus string to the phy::Nucleus state.
 * get_ion("16O6+") => Nucleus{ .A = 16, .Z = 8, .N_electrons = Some{2} }
 * get_ion("")      => Nucleus{ .A = DEFAULT_A_PRIMARY, .Z = DEFAULT_Z_PRIMARY, .N_electrons = None }
 * Second argument is only a description string for possible error message. */
::phy::Nucleus mnd::RunsheetState::get_ion(const std::string& str, const char* desc) {
	std::string_view view = mnd::trim(str);

	Option<u32> A = mnd::stou_munch(view);
	if(A.is_none()) {
		WARN("mnd::RunsheetState::get_ion(\"%s\") for input string \'%s\', "
			"couldn't parse out the mass number (A). Got `None`. "
			"Setting the value to default: %u\n", desc, str.c_str(), DEFAULT_A_PRIMARY);
		A = Some{DEFAULT_A_PRIMARY};
	}

	Option<u32> Z = phy::Z(view);
	if(Z.is_none()) {
		WARN("mnd::RunsheetState::get_ion(\"%s\") for input string \'%s\' "
			"(current view: \'%*s\'), "
			"couldn't parse out the atomic number (Z). Got `None`. "
			"Setting the value to default: %u\n", desc, str.c_str(),
			(int)view.length(), view.data(), DEFAULT_A_PRIMARY);
		Z = Some{DEFAULT_Z_PRIMARY};
	}
	
	/* Possible whitespaces before charge state */
	view = mnd::trim(view);
	
	/* Try to also parse the atomic charge.. Maybe not given. */
	Option<u32> Q = mnd::stou_munch(view);
	if(view.empty() or view[0] != '+')
		Q.reset();

	return ::phy::Nucleus {
		.A = A.unwrap(),
		.Z = Z.unwrap(),
		.N_electrons = Q.and_then([&Z](const u32 q) -> Option<u16> {
			return (Z.unwrap() >= q) ? Some{(u16)(Z.unwrap() - q)} : Some{u16{0}};
		})
	};
	/* There can be more junk after the parse, but that's fine. */
}

mnd::RunsheetState mnd::RunsheetState::from(const nlohmann::json& j) {
	return RunsheetState {
		.brho = {
			.s1_s2 = j_value<double>(j, fs::brho::s1_s2).value_or(NAN),
			.s2_s3 = j_value<double>(j, fs::brho::s2_s3).value_or(NAN),
			.s3_s4 = j_value<double>(j, fs::brho::s3_s4).value_or(NAN)
		},
		.e0 = j_value<double>(j, fs::e_primary).value_or(NAN),
		.primary  = get_ion(
			j_value<std::string>(j, fs::i_primary).value_or(""), fs::i_primary
		),
		.secondary = get_ion(
			j_value<std::string>(j, fs::i_secondary).value_or(""), fs::i_secondary
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

std::istream& phy::operator>>(std::istream& is, Nucleus& nucleus) {
	std::string token;
	if(is >> token)
		nucleus = mnd::RunsheetState::get_ion(token);

	return is;
}
