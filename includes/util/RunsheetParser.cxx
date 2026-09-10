#include "RunsheetParser.h"
#include "JSONParser.h"
#include "FromChars.h"
#include "MPhysics.hxx"
#include <string>

nlohmann::json mnd::fs::runsheet_obj {};

/* Will throw on bad parse, or if file does not exist. */
void mnd::fs::load_runsheet(const std::string& p) {
	runsheet_obj = ParseJSON(p);
}

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
	
	/* Try to also parse the atomic charge... */
	Option<u32> Q = mnd::stou_munch(view);
	if(view.empty() or view[0] != '+')
		Q.reset();

	return ::phy::Nucleus {
		.Z = Z.unwrap(),
		.A = A.unwrap(),
		.N_electrons = Q.map([](const u32 q) {
			return static_cast<uint16_t>(q);
		})
	};

	/* There can be more junk after the parse, but that's fine. */
}

mnd::RunsheetState mnd::RunsheetState::from(const nlohmann::json& j) {
	return RunsheetState {
		.brho = {
			.s1_s2 = j.value(fs::brho::s1_s2, NAN),
			.s2_s3 = j.value(fs::brho::s2_s3, NAN),
			.s3_s4 = j.value(fs::brho::s3_s4, NAN)
		},
		.e0 = j.value(fs::e_primary, NAN),
		.primary  = get_ion( j.value(fs::i_primary, ""), fs::i_primary ),
		.secondary = get_ion( j.value(fs::i_secondary, ""), fs::i_secondary )
	};
}

std::ostream& mnd::operator<<(std::ostream& os, const RunsheetState& r) {
	return os << "{ brho: {"
		<< fs::brho::s1_s2 << ": " << r.brho.s1_s2 << ", "
		<< fs::brho::s2_s3 << ": " << r.brho.s2_s3 << ", "
		<< fs::brho::s3_s4 << ": " << r.brho.s3_s4 << "}"
		<< ", {primary"
		<< ": Z: " << r.primary.Z
		<< ", A: " << r.primary.A
		<< ", Q: " << r.primary.N_electrons.map([](const auto x) { return std::to_string(x); }).value_or("none")
		<< '}'
		<< ", {secondary"
		<< ": Z: " << r.secondary.Z
		<< ", A: " << r.secondary.A
		<< ", Q: " << r.secondary.N_electrons.map([](const auto x) { return std::to_string(x); }).value_or("none")
		<< " }";
}

template<>
mnd::OptRunsheetStatePair
mnd::QueryRunsheet<
	false,
	mnd::OptRunsheetStatePair
> (std::string_view str) {
	const auto [
		basename,
		start_num,
		end_num
	] = mnd::fs::file_info(str);
	
	Option<RunsheetState> start_info = mnd::None, end_info = mnd::None;

	/* Loop over the JSON object, select the rows named as the filename
	 * base and look from there. */
	for(const auto& [_row_i, row] : fs::runsheet_obj.items()) {
		if(row.at(fs::file_name_key) != basename) continue;
				
		if(row.at(fs::start_num_key) == start_num) {
			start_info = Some{ RunsheetState::from(row) };
		}
		else if(row.at(fs::start_num_key) == end_num) {
			end_info = Some{ RunsheetState::from(row) };
		}
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
	auto [start_info, end_info] = QueryRunsheet<false, mnd::OptRunsheetStatePair>(str);

	if(start_info.is_none()) {
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename "
			"and start runnumber not found in a row in the runsheet.", (int)str.size(), str.data());
	}
	if(end_info.is_none()) {
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename "
			"and end runnumber not found in a row in the runsheet.", (int)str.size(), str.data());
	}
	if(start_info.unwrap() != end_info.unwrap()) {
		std::cerr << "mnd::QueryRunsheet: mismatch!\n"
			<< start_info.unwrap() << " ... and:\n" << end_info.unwrap() << std::endl;
		MND_THROW("mnd::QueryRunsheet(\"%*s\"): file name's inferred basename "
			"and start runnumber not found in a row in the runsheet.", (int)str.size(), str.data());
	}
	
	return start_info.unwrap();
}
