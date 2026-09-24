/* Main program: on the track multp. vs charge graph do different cuts and look
 * at the resulting spectra... */

#include "util/CLI.h"
#include "magic_enum/magic_enum.hpp"

#include "TROOT.h"
#include "TApplication.h"
#include "TParameter.h"
#include "TPaveText.h"
#include "util/Geometry.h"
#include "util/MacroHelpers.h"
#include "util/Option.hxx"
#include "util/PrettyHisto.h"
#include "util/RunsheetParser.h"
#include "util/MPhysics.h"
#include "ExperimentAssumptions.hh"
#include "util/MPhysics.h"

#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;

using namespace mnd::col;
using namespace mnd::col::literals;
using namespace phy::literals;

namespace fs = std::filesystem;
using namespace std::string_view_literals;

//#define MND_PHYSICS_PARANOIA
#define MND_MULTITHREADING_TOGGLE

using phy::AtomicNumber, phy::Nucleus, phy::RhoExpressionType;

constexpr A2 elem_to_a2(AtomicNumber e) {
	switch(e) {
		case AtomicNumber::H:  return A2{0.6, 1.5};
		case AtomicNumber::He: return A2{1.5, 2.5};
		case AtomicNumber::Li: return A2{2.5, 3.5};
		case AtomicNumber::Be: return A2{3.5, 4.5};
		case AtomicNumber::B:  return A2{4.5, 5.4};
		case AtomicNumber::C:  return A2{5.4, 6.5};
		case AtomicNumber::N:  return A2{6.5, 7.5};
		default: throw std::invalid_argument("elem_to_a2: bad atomic number?");
	}
	return A2{};
}

using rho_pad_config_t = std::tuple<
	uint32_t, // Font-type
	double,   // Font-size [inches]
	bool      // Hide statistics box
>;
using decay_text_config_t = std::tuple<
	double,  //  Font-size [inches]
	bool,    // True = do fill, No = false
	uint32_t // Fill color in hex
>;

struct SingleSelect {
	Nucleus fragment;
	u32 index;
	A2 qcut {}; // doesn't get parsed into. Is set later.
};
std::istream& operator>>(std::istream&, SingleSelect& );
std::ostream& operator<<(std::ostream&, const SingleSelect& );

template<typename T>
using Option = mnd::Option<T>;
constexpr auto None = mnd::None;

using Select = std::vector<SingleSelect>;
std::vector<phy::Nucleus> flatten_selected(const Select& );

constexpr mnd::col::RGBA DEFAULT_FILL_COL_EX = 0x91E965_c + 0.77_o;

int main(int argc, char* argv[]) {
	CLI::App app{"This program analyses the Hit level file and produces (hopefully) decent physics."};

	std::vector<std::string> fileName, info;
	double distance_cut = NAN;
	Option<A2> sci21_cut;
	Option<A2> sci22_cut;
	Option<A2> sci31_cut;
	size_t max_events = -1;
	u32 nthreads = 1;

	Select selected{};
	bool do_scaling = false;
	bool fmt_isotope = false;
	std::vector<canvas::Extension> save = {};

	A3 theta_binning = {200,0,200};
	A3 rho_binning = {300,0,300};
	
	rho_pad_config_t rho_config{ 132, 0.045, false };
	/* Few params really only related to the pyplot histogram. */
	std::string rho_title, rho_xlabel;
	Option<double> py_linewidth;
	Option<RGBA> py_fillcol = mnd::Some{ DEFAULT_FILL_COL_EX };
	Option<RGBA> py_linecol;

	add_logged_option(app, "-f,--file", fileName, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	add_logged_option<DisplayDefault::No>(app, "-m,--max-events", max_events,
		"Max events (at most) taken from each ROOT file. Default: all entries.");
	add_logged_option(app, "-s,--select", selected, "Select which fragments to explicity gate upon, delimited by ':'")
		->delimiter(':')
		->type_name("ELEM,N | A,B,N");
	add_logged_option<DisplayDefault::No>(app, "-d,--distance", distance_cut,
		"Distance from vertex cut (all tracks must be below this threshold). Default no cut. "
		"NOTE: not assigning this field will default to simple angles between all FOOT tracks.")
		->check(CLI::Range(0.0, 2.0));
	add_logged_option(app, "--sci21", sci21_cut, "SCI21 QDC cut (also implying multiplicity 1).")
		->delimiter(',');
	add_logged_option(app, "--sci22", sci22_cut, "SCI22 QDC cut (also implying multiplicity 1).")
		->delimiter(',');
	add_logged_option(app, "--sci31", sci31_cut, "SCI31 QDC cut (also implying multiplicity 1).")
		->delimiter(',');
	add_logged_flag(app, "--scaling", do_scaling, "Rescale the ρ angle based on the "
		"calculated kinetic energy of the reaction inside the target");
	add_logged_option(app, "-o,--save", save, "Save the resulting canvases as one or more extensions.")
		->delimiter(',');
	add_logged_option(app, "-i,--info", info,
		"Sequence of strings, such that the output will be saved descending these directories in autosave.")
		->delimiter(',');
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
#ifndef MND_MULTITHREADING_TOGGLE
		BOLD " Parameter ignored. Compiled purely singlethreaded." KNRM
#endif // !MND_MUL
		)->check(CLI::PositiveNumber);
	add_logged_option(app, "--rho-binning", rho_binning, "Binning of the ρ-value, in mrad.")
		->delimiter(',');
	add_logged_option(app, "--theta-binning", theta_binning,
		"Binning of the θ values, in mrad. This represents the angle between heavy fragment and all the lighter ones.")
		->delimiter(',');
	add_logged_option(app, "--rho-title", rho_title, "Title of the ρ-value histogram. Can be latex'ed (inside \'\' block)");
	add_logged_option(app, "--rho-xlabel", rho_xlabel, "X-axis label of the ρ-value histogram. Can be latex'ed (inside \'\' block)");
	add_logged_option(app, "--rho-tstyle", rho_config, "Configure the ROOT style of ρ-value 1D histogram.\n"
		"Arg[0]: font style.\nArg[1]: font width (in inches).\nArg[2]: true => show statistics box.")
		->delimiter(',');
	add_logged_option(app, "--py-linewidth", py_linewidth, "Linewidth for the ρ-value histogram to be exported from Python.");
	add_logged_option(app, "--py-fillcol", py_fillcol, "Line col (ARGB) for the ρ-value histogram to be exported from Python. "
		"By default, taken from original TH1P")
		->default_str( mnd::to_string(DEFAULT_FILL_COL_EX) );
	add_logged_option(app, "--py-linecol", py_linecol, "Line col (ARGB) for the ρ-value histogram to be exported from Python. "
		"If left as none, taken from original TH1P");
	add_logged_flag(app, "--fmt-isotope", fmt_isotope,
		"Add Isotope indicator inside the autogenerated LaTeX labels for ρ-value");

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;
	if(fileName.size() == 0)
		ERROR("To continue, must supply at least one file name!\n");

	mnd::fs::load_runsheet();
	std::unique_ptr<mnd::RunsheetState> run_info;
	{
		auto sFront = std::make_unique<mnd::RunsheetState>( mnd::QueryRunsheet(fileName.front()) );
		auto sBack  = std::make_unique<mnd::RunsheetState>( mnd::QueryRunsheet(fileName.back()) );
		/* I "trust" the OS that the allocation won't fail, e.g. I'm not dereffing a null here :-) */
		if(*sFront != *sBack) {
			WARN("Runsheet status for initial file '%s' is: \n", fileName.front().c_str());
			std::cerr << *sFront << std::endl;
			fprintf(stderr, ".. and for the final file '%s' is: \n", fileName.back().c_str());
			std::cerr << *sBack << "\n" KBH_RED ">> does not match!\n" KNRM;
		}
		run_info = std::move(sFront);
	}
	WARN("Run-info parsed as: "); std::cerr << *run_info << std::endl;
	WARN("Primary beam: "); std::cerr << run_info->primary << std::endl;
	WARN("Secondary beam: "); std::cerr << run_info->secondary << std::endl;
	if(!run_info->primary.valid() or !run_info->secondary.valid())
		ERROR("Either primary or secondary beam is invalid. Parse from runsheet was unsuccessful.\n");
	
	const double A_beam = run_info->secondary.A;
	const double Z_beam = run_info->secondary.Z;

	/* Average kinetic energy per nucleon, as beam enters the S2. */
	const double ekin_s2 = phy::EKin(
		A_beam, Z_beam,
		phy::Brho_t{run_info->brho.s1_s2}
	);
	/* Average kinetic energy per nucleon, just before the target. */
	const double avg_ekin_before_target =
		ekin_s2 - mnd::assume::s2::loss_upto_target;
	/* Average kinetic energy per nucleon, just after the target. */
	const double avg_ekin_after_target =
		ekin_s2 - mnd::assume::s2::loss_upto_target - mnd::assume::s2::loss_in_target;

	/* Average kinetic energy at the centre of the target. */
	const double avg_ekin_reaction = (avg_ekin_before_target + avg_ekin_after_target) / 2.0;
	/* Average beta corresponding to this energy. */
	const double beta_nominal = phy::Beta(
		A_beam, Z_beam,
		phy::EKin_t{avg_ekin_reaction}
	);
	const double gamma_nominal = phy::Gamma(beta_nominal);
	const double beta_gamma_n = beta_nominal * gamma_nominal;

	WARN("Assumed values for just before/after 9Be target: "
		MND_RGB_COL(250, 250,  70) "%.2f MeV/u" KNRM
		" and "
		MND_RGB_COL( 70, 160, 250) "%.2f MeV/u\n" KNRM,
		avg_ekin_before_target, avg_ekin_after_target);

	mnd::python::poke();
	ROOT::EnableThreadSafety();

	TApplication rootApp("app", 0, 0);

	using namespace mnd::geom;

	double Cr, Cq, Ct, max_cost, max_cost_f;
	double target_depth = 0.0;
	{
		const auto& fname = fileName.front();
		std::array<double, 3>* c;
		TParameter<double>* m;
		auto f = std::make_unique<TFile>(fname.c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2);
		max_cost = m->GetVal();
		get_obj(f, m, "FOOT_max_cost_f");
		max_cost_f = m->GetVal();
		FOOTBoxParam* box;
		get_obj(f, box, "FOOT_box");
		target_depth = box->target.Width();
	}
	constexpr auto TARGET_Z = TFOOTHitCont::TARGET_Z;
	WARN("9Be S2 target's outer edge placed nominally at " MND_RGB_COL(255,128,0)
		"z₀ = %.2f mm" KNRM " with width " MND_RGB_COL(255,128,0)
		"w = %.2f mm" KNRM"\n", TARGET_Z, target_depth);
	
	const double loss_rate_per_depth =
		(avg_ekin_after_target - avg_ekin_before_target) / target_depth;

	if(do_scaling) {
		if(target_depth <= 0.0)
			ERROR("Scaling requested (-s,--scaling) but target width "
				"calculated as %.2f <= 0.0 ? Not allowed!\n", target_depth);

		WARN("Normalising will slightly scale ρ-value, per entry, to the "
			"average energy in the center of the 9Be target: "
			KBH_CYN "%.3f" KNRM " , which means: "
			KBH_BLU "beta = %.5f. " KNRM "\nEnergy loss rate per depth calculated as: "
			MND_RGB_COL(255,128,0) "%.2f AMeV/mm" KNRM " while per x distance: "
			MND_RGB_COL(128,255,0) "%.3f AMeV/mm" KNRM "\n",
			avg_ekin_reaction, beta_nominal, loss_rate_per_depth, mnd::assume::s2::e_dispersion[1]);
	}

	/* Sanitize some CLI passed in arguments... */
	const u32 sum_n_tracks_required = std::accumulate(
		selected.begin(), selected.end(), 0, [](u32 sum, const auto& x) {
			return sum + x.index;
		}
	);
	if(sum_n_tracks_required == 1)
		ERROR("Cannot require just a single track to form a valid vertex. -s,--select option "
		      "must require at least two particles.\n");

	/* Populate the missing field of the SingleSelect entry with the charge range. */
	for(auto& [nucleus, n_required, qcut] : selected) {
		if(n_required == 0)
			ERROR("Must require 1 or more tracks, cannot set to require 0 tracks. Stop that.");
		qcut = elem_to_a2( nucleus.ChemElem().unwrap() );
	}
	
	/* Sort out the selected nuclei, in descending charge.
	 * Convention from RNFOOTTrack::operator<(...)
	 * If charges are the same, then compare the mass.
	 * Also use this time to do some sanity checks for the ranges. */
	std::sort(selected.begin(), selected.end(),
		/* Sort in descending order, according to the lower range boundary. */
		[](const SingleSelect& lhs, const SingleSelect& rhs) {
			if(lhs.qcut[0] > lhs.qcut[1]) {
				WARN("Range: "); std::cerr << lhs.fragment << " qcut range isn't valid (lhs >= rhs)? Is: " << lhs.qcut;
				throw std::runtime_error("Err");
			}
			if(rhs.qcut[0] > rhs.qcut[1]) {
				WARN("Range: "); std::cerr << rhs.fragment << " qcut range isn't valid (rhs >= rhs)? Is: " << rhs.qcut;
				throw std::runtime_error("Err");
			}
			return (lhs.fragment.Z == rhs.fragment.Z)
			     ? (lhs.fragment.mass() > rhs.fragment.mass())
			     : (lhs.fragment.Z > rhs.fragment.Z);
		}
	);
	if( std::adjacent_find(selected.begin(), selected.end(),
		[](const SingleSelect& lhs, const SingleSelect& rhs) {
			return lhs.qcut[0] < rhs.qcut[1];
			// If LHS elem's lower bound is smaller than RHS elem's upper bound, means their intervals overlap.
		}
	) != selected.end() ) {
		ERROR("Two intervals from the 'selected' (-s,--select) field are overlapping. Not allowed");
	}
	
	RhoExpressionType rho_ex_type = RhoExpressionType::full;
	/* Sanitize the angle part. Check if either the reaction of interest is:
	 * [1] X => Y + Np
	 * [2] X => y_1 + y_2 + ... + y_N 
	 * [3] X => N y
	 * [4] UNKNOWN , in case no selection (-s,--select) was given. */
	if(selected.size() == 2 &&
	   selected.front().index == 1 &&
	  (selected.back().fragment == "1H"_n ||
	   selected.back().fragment == "1H1+"_n)) {
		rho_ex_type = RhoExpressionType::p;
	}
	else if(selected.size() == 1) {
		rho_ex_type = RhoExpressionType::equinuclear;
	}
	else if(selected.empty()) {
		rho_ex_type = RhoExpressionType::unknown;
	}

	/* Sequence heavier-to-lighter of daughter nuclei. */
	const std::vector<Nucleus> nuclei = flatten_selected(selected);
	
	/* From the selection, calculate the compound nucleus. */
	Nucleus mother{.A = 0, .Z = 0, .N_electrons = None };
	for(const auto& nuc: nuclei) {
		mother &= nuc;
	}
	{
		std::string_view rho_type_label =  magic_enum::enum_name(rho_ex_type);
		WARN("Matched the ρ-analysis type: it is: %s%.*s%s, the decay of "
			MND_RGB_COL(141,249,155) "%s" KNRM "\n",
			BOLD, (int)rho_type_label.length(), rho_type_label.data(), KNRM,
			mother.to_string().c_str());
	}

	if(mother.valid() and nuclei.size() > 1) {
		WARN(BOLD "The complete requested reaction deduced as: "
			KNRM "%s" MND_RGB_COL(255,65,224) "\n",
			phy::format_reaction(phy::Nucleus::Represent::Normal,
				std::vector{mother}, nuclei ).c_str()
		);
		WARN("Daughter nuclei: "); std::cerr << nuclei << std::endl;
	}
	else if(!mother.valid() and nuclei.size() > 1) {
		ERROR("Requesting a reaction channel for mother nucleus: %s , but "
			"the nucleus qualified as invalid. Maybe its mass not given?\n", 
			mother.to_string().c_str());
	}
	
	/* Why the f*** ROOT can't parse a LaTeX expression naturally? Amazing job, devs. :b */
	std::string rho_expression_rootex = "";
	std::string rho_expression_pythex = "$";
	switch(rho_ex_type) {
		case(RhoExpressionType::p): {
			const auto& [heavy, nheavy_, _] = selected.front();
			const auto& [proton, nprotons, __] = selected.back();
			if(nheavy_ != 1)
				ERROR("Something is wrong. Heavy ion: '%s' selected as quantity %u and not 1\n?",
					heavy.to_string().c_str(), nheavy_);
			
			const auto heavy_ion_label_rootex = heavy.chem_to_string(phy::Nucleus::Represent::Rootex, fmt_isotope);
			const auto heavy_ion_label_pythex = heavy.chem_to_string(phy::Nucleus::Represent::Pythex, fmt_isotope);
			if(nprotons > 1) {
				rho_expression_rootex += mnd::msg("#frac{1}{#sqrt{1 + %u m[p] / m[%s]}}"
					" #sqrt{", nprotons, heavy_ion_label_rootex.c_str());
				rho_expression_pythex += mnd::msg("\\frac{1}{ \\sqrt{1 + %u \\frac{ m_{\\mathrm{p}} }{ m_{%s} } } }"
					" \\sqrt{", nprotons, heavy_ion_label_pythex.c_str());

				for(u32 i=1; i <= nprotons; ++i) {
					rho_expression_rootex += Form("%s#theta_{%s-p_{%u}}^{2} ", (i==1)? "": "+ ",
						heavy_ion_label_rootex.c_str(), i);
					rho_expression_pythex += Form("%s\\theta_{%s-\\mathrm{p}_%u}^2", (i==1)? "": "+",
						heavy_ion_label_pythex.c_str(), i);
				}
				rho_expression_rootex += mnd::msg("+ #frac{m[p]}{m[%s]} ", heavy_ion_label_rootex.c_str());
				rho_expression_pythex += mnd::msg("+ \\left(\\frac{m_\\mathrm{p}}{m_{%s}}\\right)", heavy_ion_label_pythex.c_str());
				if(nprotons > 2) {
					rho_expression_rootex += '(';
					rho_expression_pythex += "\\left(";

				}
				for(u32 i=1; i <= nprotons; ++i) {
					for(u32 j=i+1; j <= nprotons; ++j) {
						rho_expression_rootex += mnd::msg("%s#theta_{p_{%u} - p_{%u}}^{2} ",
							(i==1 && j==2)? "": "+ ", i, j);
						rho_expression_pythex += mnd::msg("%s\\theta_{\\mathrm{p}_{%u} - \\mathrm{p}_{%u}}^2",
							(i==1 && j==2)? "": "+", i, j);
					}
				}
				if(nprotons > 2) {
					rho_expression_pythex += "\\big)";
					rho_expression_rootex += ")";
				}
				rho_expression_pythex += '}';
				rho_expression_rootex += '}';
			}
			else {
				rho_expression_rootex += mnd::msg("#frac{1}{#sqrt{1 + m[p]/m[%s]}}"
					" #theta_{%s-p}",
					heavy_ion_label_rootex.c_str(), heavy_ion_label_rootex.c_str());
				rho_expression_pythex += mnd::msg("\\frac{1}{ \\sqrt{1 + \\frac{ m_{\\mathrm{p}} }{ m_{%s} } } }"
					"\\,\\theta_{%s-\\mathrm{p}}",
					heavy_ion_label_pythex.c_str(), heavy_ion_label_pythex.c_str());
			}
			break;
		}
		case(RhoExpressionType::equinuclear): {
			const auto& [nuc, nitems, _] = selected.back();
			const auto nuc_label_rootex = nuc.chem_to_string(phy::Nucleus::Represent::Rootex, fmt_isotope);
			const auto nuc_label_pythex = nuc.chem_to_string(phy::Nucleus::Represent::Pythex, fmt_isotope);

			rho_expression_rootex += mnd::msg("#sqrt{#frac{m[%s]}{%u m[p]} #sqrt{",
				nuc_label_rootex.c_str(), nitems);
			rho_expression_pythex += mnd::msg("\\sqrt{\\frac{m_{%s}}{%u m_{\\mathrm{p}} }} \\sqrt{",
				nuc_label_pythex.c_str(), nitems);

			for(u32 i=1; i<=nitems; ++i) {
				for(u32 j=i+1; j<=nitems; ++j) {
					rho_expression_rootex +=  mnd::msg("%s#theta_{%s(%u) - %s(%u)}^{2}",
						(i==1 && j==2)? "": "+",
						nuc_label_rootex.c_str(), i, nuc_label_rootex.c_str(), j);
					rho_expression_pythex +=  mnd::msg("%s\\theta_{%s(%u) - %s(%u)}^2",
						(i==1 && j==2)? "": "+",
						nuc_label_pythex.c_str(), i, nuc_label_pythex.c_str(), j);
				}
			}
			rho_expression_pythex += '}';
			rho_expression_rootex += '}';
			break;
		}
		case(RhoExpressionType::full): {
			const auto mother_label_rootex = mother.chem_to_string(phy::Nucleus::Represent::Rootex, fmt_isotope);
			const auto mother_label_pythex = mother.chem_to_string(phy::Nucleus::Represent::Pythex, fmt_isotope);
			rho_expression_rootex += mnd::msg(
				"#sqrt{#sum_{i=1}^{%u} #sum_{j=i+1}^{%u} #frac{m_{i} m_{j}{ M[%s] m[p] } #theta_{ij}^{2} }",
					sum_n_tracks_required, sum_n_tracks_required, mother_label_rootex.c_str());
			rho_expression_pythex += mnd::msg(
				"\\sqrt{\\sum_{i=1}^{%u} \\sum_{j=i+1}^{%u} \\frac{m_i m_j}{ M_{%s} m_{\\mathrm{p}} } \\theta_{ij}^2 }",
					sum_n_tracks_required, sum_n_tracks_required, mother_label_pythex.c_str());
			break;
		}
		case(RhoExpressionType::unknown): {
			rho_expression_rootex += "#rho^{2} #approx #sum_{i=1}^{N} #theta_{HI, i}^{2}";
			rho_expression_pythex += "\\rho \\simeq \\sqrt{\\sum_{i=1}^{N} \\theta_{\\mathrm{HI}, \\mathrm{p}}^2}";
			break;
		}
	}
	rho_expression_pythex += '$';

	WARN("Rho expression in ROOT: " BOLD "%s" KNRM "\n", rho_expression_rootex.c_str());
	WARN("Rho expression in Pyth: " BOLD "%s" KNRM "\n", rho_expression_pythex.c_str());
	if(rho_title.empty())
		rho_title = rho_expression_pythex;

	const bool vertex_dist_cut_given = std::isfinite(distance_cut);

	auto h1_track_mult = TH1P{"Track multiplicity [unique tracks]@-1 means no FOOT in event.", kRed-1, 11, -1.5, 9.5};
	if(vertex_dist_cut_given and sum_n_tracks_required > 0)
		h1_track_mult.AppendToTitle(" 0 means event not viable.");

	auto h2_q_vs_mult = TH2P{"Track charge [Q]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 40, 0,8};
	auto h2_score_vs_mult = TH2P{"Track score [a.u.]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 300, 0, 50};
	auto h2_track_distance = TH2P{"Track distances to vertex [mm]:Track I@smaller x-axis means larger Q particle",
		(vertex_dist_cut_given && sum_n_tracks_required>0)? sum_n_tracks_required: 10,
		-0.5,
		(vertex_dist_cut_given && sum_n_tracks_required>0)? sum_n_tracks_required-0.5: 9.5,
		600, 0, vertex_dist_cut_given? (2*distance_cut): 30.0};

	auto h1_track_angle = TH1P{"#theta(HI,LI) [mrad]@Between one heavy and all lighter tracks", kMagenta+1, 200, 0, 100};
	auto h1_angle_ex = TH1P{"((h1_angle_ex)) #rho [mrad]@#rho value", DEFAULT_FILL_COL_EX,
		rho_binning[0], rho_binning[1], rho_binning[2]
	};
	h1_angle_ex->SetTitle( rho_expression_rootex.c_str() );
	auto h2_vertex_z = TH2P{"#rho [mrad]:Vertex z [mm]@Traced by the FOOT", 140, -80, 80,
		rho_binning[0], rho_binning[1], rho_binning[2]};
	auto h2_vertex_x = TH2P{"#rho [mrad]:Vertex x [mm]@Traced by the FOOT", 140, -60, 60,
		rho_binning[0], rho_binning[1], rho_binning[2]};
	if(do_scaling) {
		h1_angle_ex.AppendToAxisTitle(" (beta corrected)"sv);
		h2_vertex_z.AppendToAxisTitle(" (beta corrected)"sv);
		h2_vertex_x.AppendToAxisTitle(" (beta corrected)"sv);
	}
	auto h2_nlayers_hit = TH2P{"N layers hit in a track:Track Ip@Full FOOT system, smaller x-axis means larger Q particle",
		10, -0.5, 9.5,   RNFOOTHit::N_PAIRS, 0.5, RNFOOTHit::N_PAIRS+0.5 };
	auto h2_rho_vs_theta = TH2P{"#rho [mrad]:#theta(p,HI) [mrad]@Traced by the FOOT",
		theta_binning[0], theta_binning[1], theta_binning[2],
		rho_binning[0], rho_binning[1], rho_binning[2]};
	if(rho_ex_type ==  RhoExpressionType::p) {
		const Nucleus& hi = nuclei.front();
		h2_rho_vs_theta->GetXaxis()->SetTitle(
			Form("#theta(p, %s)", hi.chem_to_string(phy::Nucleus::Represent::Rootex, fmt_isotope).c_str())
		);
	}
	auto h1_sci21 = TH1P{"SCI21 QDC mean [QDC units]", 0xCB00CB_c, 500, 300, 4000};
	auto h1_sci22 = TH1P{"SCI22 QDC mean [QDC units]", 0x0070DD_c, 500, 300, 4000};
	auto h1_sci31 = TH1P{"SCI31 QDC mean [QDC units]", 0x009B2F_c, 500, 300, 4000};
	auto h1_sci21_cut  = TH1P{"((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000};
	auto h1_sci22_cut  = TH1P{"((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000};
	auto h1_sci31_cut  = TH1P{"((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, 500, 300, 4000};
	auto h1_sci21_cut2 = TH1P{"((h1_cut2)) SCI21 QDC mean [QDC units]@With cut and FOOT selection", 0x890389_c, 500, 300, 4000};
	auto h1_sci22_cut2 = TH1P{"((h1_cut2)) SCI22 QDC mean [QDC units]@With cut and FOOT selection", 0x6180FD_c, 500, 300, 4000};
	auto h1_sci31_cut2 = TH1P{"((h1_cut2)) SCI31 QDC mean [QDC units]@With cut and FOOT selection", 0x7DE69D_c, 500, 300, 4000};

	show_console_cursor(false);

#ifdef MND_MULTITHREADING_TOGGLE
	DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileName, nthreads, [=, &progress](size_t i, auto fname) mutable {
#else
	for(size_t i=0; i<fileName.size(); ++i) {
		const auto& fname = fileName[i];
#endif
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTHit>("FOOT");
		auto frs = model->MakeField<RNFRSHit>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);
		const size_t nentries = ( (max_events < ntuple->GetNEntries()) ? max_events : ntuple->GetNEntries() );

#ifdef MND_MULTITHREADING_TOGGLE
		auto _bar_ptr = std::make_unique<ProgressBar>
#else
	auto bar = ProgressBar
#endif
		(
			option::BarWidth{55},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("%zu/%zu: %'zu (%s)", i+1, fileName.size(), nentries, fname.c_str())},
			option::ForegroundColor{ indicators::next_col() },
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		);

#ifdef MND_MULTITHREADING_TOGGLE
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];
#endif
	
		for(size_t entryId{0}; entryId < nentries; ++entryId ) {
#ifdef MND_MULTITHREADING_TOGGLE
			if( mnd::PrintProgress(bar, entryId, nentries, 1000) )
				progress.print_progress();
#else
			mnd::PrintProgress(bar, entryId, nentries, 1000);
#endif

			ntuple->LoadEntry(entryId);

			const auto& sci21 = frs->cal.sci[0];
			const auto& sci22 = frs->cal.sci[1];
			const auto& sci31 = frs->cal.sci[2];
			
			if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);
			
			/* Skip the event entirely in case a SCI cut isn't met. */
			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() >= 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31_cut->Fill(sci31.E);
			
			bool valid_vertex_found = false;

			/* In this case, don't cut on any charges etc,.. just take the whole event and try to do
			* general vertexing, angles, some kind of rho, etc. */
			if(!vertex_dist_cut_given) {
				const size_t N = foot->track.size();
				h1_track_mult->Fill( foot->HasData()? N: -1 );

				std::vector<Line3D> tracks {};

				/* If no selection is provided, then just fetch everything. */
				for(size_t i=0; i<N; ++i) {
					const RNFOOTTrack& t = foot->track[i];

					h2_q_vs_mult->Fill(N, t.Q);
					h2_score_vs_mult->Fill(N, t.score);

					tracks.push_back(*t);
				}
				if(N < 2) continue;

				const Point3D vertex = FindVertex(tracks);
				for(size_t i=0; i<N; ++i) {
					h2_track_distance->Fill( i, tracks[i].DistanceTo(vertex) );
				}
				double sum2 = 0;
				double theta;
				const Line3D heavy_track = tracks.front();
				for(size_t i=1; i < N; ++i) {
					theta = phy::MRAD_CVT * tracks[i].AngleRelativeTo( heavy_track );
					sum2 += theta*theta;
					h1_track_angle->Fill( theta );
				}
				const double invariant_theta = sqrt(sum2);
				h1_angle_ex->Fill( invariant_theta );
				h2_vertex_z->Fill( vertex.z, invariant_theta);
				h2_vertex_x->Fill( vertex.x, invariant_theta);

				valid_vertex_found = !vertex.is_null();
			}
			else { while(foot->track.size() >= 2) {
				VertexingResult<RNFOOTTrack> vtr =
					FindVertexingTracksMut(foot->track, distance_cut);

				/* Can promptly break from the `while`; no further tracks can make a vertex anymore. */
				if(!vtr.valid()) break;

				const std::vector<RNFOOTTrack>& ftracks = vtr.tracks;
				const Point3D& vertex = vtr.vertex;
				const u32 n_tracks_selected = vtr.tracks.size();
				
				if(sum_n_tracks_required > 0 and (n_tracks_selected != sum_n_tracks_required))
					continue;
				
				u32 mask = (1U << n_tracks_selected) - 1; /* sequence 0b 1111.. */
				u32 unique_track_index = 0;
				/* Next, we require that this vertexing sequence of tracks also
				* satisfies the charge cut. Selection windows are already
				* presorted in descending charge. */
				for(const auto& [nuc, n_tracks_required, charge_cut] : selected) {
					for(u32 i=0; i<n_tracks_required; ++i) {
						if(mnd::IsInside(
							ftracks.at(unique_track_index).Q,
							charge_cut)
						) { // Reset the unique bit.
							mask &= ~(1U << unique_track_index);
						}
						++unique_track_index;
					}
				}
				#ifdef MND_PHYSICS_PARANOIA
				assert( unique_track_index == sum_n_tracks_required
					&& "Huh seriously am I so bad?" );
				#endif

				/* `mask` now must be 0, otherwise we didn't catch all unique tracks for this
				* specific charge interval. */
				if(sum_n_tracks_required > 0 && mask != 0)
					continue;
			
				/* Due to lexicographical combinations, the cute fact is that
				* vtr.tracks is also already sorted accordingly :) */
				#ifdef MND_PHYSICS_PARANOIA
				assert( std::is_sorted(ftracks.begin(), ftracks.end())
					&& "Huh really?" );
				#endif
				
				/* At this point, selection is completed.
				* Manifest back the Line3D objects to do ρ-value calculations. */
				std::vector<Line3D> tracks;
				tracks.reserve(n_tracks_selected);
				for(const auto& ft : ftracks) {
					h2_q_vs_mult->Fill(n_tracks_selected, ft.Q);
					h2_score_vs_mult->Fill(n_tracks_selected, ft.score);
					tracks.push_back(*ft);
				}
				
				double rho_val = phy::rho(
					rho_ex_type,
					nuclei,
					tracks
				);
				if(do_scaling) {
					double depth = std::clamp(
						vertex.z,
						TARGET_Z - target_depth,
						TARGET_Z
					) - (TARGET_Z - target_depth);
					double x_ = mnd::clamp(
						vertex.x,
						mnd::assume::s2::x_bound
					);
					const double ekin_per_n = avg_ekin_before_target
						+ loss_rate_per_depth * depth
						+ mnd::assume::s2::e_dispersion[0]
						+ mnd::assume::s2::e_dispersion[1] * x_;

					const double beta_gamma = phy::BetaGamma(
						A_beam, Z_beam, // secondary beam runinfo asserted to be valid
						phy::EKin_t{ekin_per_n}
					);
					rho_val *= (beta_gamma / beta_gamma_n);
				}
				
				h1_angle_ex->Fill( rho_val );
				h2_vertex_z->Fill( vertex.z, rho_val );
				h2_vertex_x->Fill( vertex.x, rho_val );

				for(u32 i=0; i<n_tracks_selected; ++i) {
					h2_track_distance->Fill( i, tracks[i].DistanceTo(vertex) );
					h2_nlayers_hit->Fill( i, ftracks[i].n );
				}

				const mnd::geom::Line3D& heavy_track = tracks.front();
				for(u32 i=1; i < n_tracks_selected; ++i) {
					double theta = phy::MRAD_CVT * tracks[i].AngleRelativeTo( heavy_track );
					h1_track_angle->Fill( theta );
					h2_rho_vs_theta->Fill(theta, rho_val);
				}
				h1_track_mult->Fill(n_tracks_selected);
				valid_vertex_found = true;

			} /* while(...) */ } // if(vertex_dist_cut_given)
			
			if(valid_vertex_found and vertex_dist_cut_given) {
				h1_sci21_cut2->Fill(sci21.E);
				h1_sci22_cut2->Fill(sci22.E);
				h1_sci31_cut2.Fill(sci31.E);
			}
		} // for(size_t entryId{0}; entryId < nentries; ++entryId )

#ifdef MND_MULTITHREADING_TOGGLE
		bar.mark_as_completed();
#endif
	}

#ifdef MND_MULTITHREADING_TOGGLE
); // parallel_process
#endif

	show_console_cursor(true);

	TCanvas* cm = new TCanvas("Multp", "Different tracks and extra info", 2150, 1400);
	cm->Divide(3,2);
	cm->cd(1); h2_q_vs_mult.Draw("COLZ"); gPad->SetLogz();
	cm->cd(2); h2_score_vs_mult.Draw("COLZ"); gPad->SetLogz();
	cm->cd(3); h2_track_distance.Draw();
	cm->cd(4); h1_track_mult.Draw();
	cm->cd(5); h2_nlayers_hit.Draw("COLZ");
	cm->cd(6); new PLatex(0.06,
		"Coefficients: ",
		Form("Cr = %.1f mm^-2", Cr),
		Form("Cq = %.1f Q^-2", Cq),
		Form("Ct = %.1f mm^-2", Ct),
		Form("max cost for candidate: %.1f", max_cost),
		Form("max cost for whole track: %.1f", max_cost_f)
	);
	
	/* For the ρ-value histogram, add upper axis that is reflective on the transversal
	 * part of the decay kinetic energy:
	 *   ρ(T⟂) = 1/(beta*gamma) sqrt(2*T⟂ / (mp*c^2))
	 *   T⟂(ρ) = mp*c^2 / 2 * (beta*gamma)^2 * ρ^2
	 */
	const double C_Tl_to_rho = /* Units [rad / sqrt(MeV)] */
		1.0 / beta_gamma_n * sqrt(2.0 / (phy::mp * phy::nuc::c*phy::nuc::c ));
	const double C_rho_to_Tl = /* Units [MeV / rad^2] */
		(phy::mp * phy::nuc::c*phy::nuc::c )/2 * (beta_gamma_n * beta_gamma_n);

	auto f_rho_to_Tl =
	[C = C_rho_to_Tl](double rho /* [mrad] */) -> double /* [MeV] */ {
		double rho_rad = rho / phy::MRAD_CVT;
		return C * rho_rad * rho_rad;
	};
	auto f_Tl_to_rho =
	[C = C_Tl_to_rho](double Tl /* [MeV] */) -> double /* [mrad] */ {
		return phy::MRAD_CVT * C * std::sqrt(Tl);
	};
	TCanvas* ct = new TCanvas("Physics", "Recognized tracks and angles", 2150, 1400);
	ct->Divide(2,2);
	ct->cd(1); h2_rho_vs_theta.Draw("COLZ");
	ct->cd(2); h2_vertex_x.Draw<TH2P::Y>(
		f_rho_to_Tl,
		f_Tl_to_rho,
		"Q_{#perp}  [MeV]",
		"COLZ"
	);
	ct->cd(3); h1_angle_ex.Draw(
		f_rho_to_Tl,
		f_Tl_to_rho,
		"Q_{#perp}  [MeV]",
		"HIST"
	);
	h1_angle_ex->SetStats(std::get<2>(rho_config));
	if(auto* title = dynamic_cast<TPaveText*>(gPad->GetPrimitive("title"))) {
		title->SetTextFont( std::get<0>(rho_config) );
		title->SetTextSize( std::get<1>(rho_config) );
		title->SetBorderSize(6);
		title->SetShadowColor( (0x343434_c).GetColorCode() );
	}
	ct->cd(4); h2_vertex_z.Draw<TH2P::Y>(
		f_rho_to_Tl,
		f_Tl_to_rho,
		"Q_{#perp}  [MeV]",
		"COLZ"
	);

	TCanvas* ct2 = new TCanvas("Physics2", "Recognized tracks and more", 2150, 1400);
	ct2->Divide(2,2);
	ct2->cd(1); h1_track_angle.Draw();
	ct2->cd(2); h2_track_distance.Draw("COLZ");
	/* TODO: 2 more slots here.. */

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2150, 1400);
	cs->Divide(3,3);
	cs->cd(1); h1_sci21.Draw();
	cs->cd(2); h1_sci22.Draw();
	cs->cd(3); h1_sci31.Draw();
	cs->cd(4); h1_sci21_cut.Draw();
	cs->cd(5); h1_sci22_cut.Draw();
	cs->cd(6); h1_sci31_cut.Draw();
	cs->cd(7); h1_sci21_cut2.Draw();
	cs->cd(8); h1_sci22_cut2.Draw();
	cs->cd(9); h1_sci31_cut2.Draw();

	WARN("Info: "); std::cerr << MND_RGB_COL(140,85,255) << info << KNRM "\n";

	canvas::save_all<canvas::Exe>( save, mnd::to_views(info) );

	auto py_histstyle = mnd::plot::HistStyle{}.stairs();
	if(py_fillcol.is_some()) {
		py_histstyle = std::move(py_histstyle)
			.fill()
			.facecolor(py_fillcol.unwrap());
	}
	if(py_linewidth.is_some()) {
		py_histstyle = std::move(py_histstyle)
			.line_width(py_linewidth.unwrap());
	}
	if(py_linecol.is_some()) {
		py_histstyle = std::move(py_histstyle)
			.edgecolor(py_linecol.unwrap());
	}

	mnd::plot::Figure {}
		.plot(*h1_angle_ex, py_histstyle)
		.save_dpi(200)
		.xlabel(
			 !rho_xlabel.empty()
			? rho_xlabel
			: R"($\rho\,[\mathrm{mrad}]$)")
		.ylabel(&h1_angle_ex.h)
		.grid()
		.title(
			 !rho_title.empty()
			? rho_title
			: R"($\rho = \sqrt{\sum_{i} \theta_{\mathrm{HI}-p_{i}^2}$)")
		.save(
			fs::path{"autosave"}
			/ mnd::fs::current_executable_name()
			/ "py"
			/ mnd::fs::path_sequence(info)
			/ "excitation.png"
		);

	WARN("End-of-main\n");
	rootApp.Run(); return 0;
}

/* Target format: "4He,2" */
std::istream& operator>>(std::istream& is, SingleSelect& out) {
	std::string s;
	is >> s;
	auto parts = mnd::split_view(s, ','); // vector<string_view>

	switch(parts.size()) { // {elem, n}, e.g.: { "4He"sv, "2"sv }
		case(2) : {
			Nucleus nuc;
			try {
				nuc = phy::Nucleus::get_ion(std::string{parts[0]}, true, "std::istream& operator>>");
			} catch(std::exception& e) {
				YELL("%s\n", e.what());
				is.setstate(std::ios::failbit);
				return is;
			}

			Option<u32> n = mnd::stou(parts[1]);
			if(n.is_none()) {
				is.setstate(std::ios::failbit);
				return is;
			}

			out = SingleSelect {
				.fragment = std::move(nuc),
				.index = n.unwrap()
			};
			break;
		}
		default:
			is.setstate(std::ios::failbit);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, const SingleSelect& x) {
	return os << x.fragment << ':' << x.index;
}

std::vector<phy::Nucleus> flatten_selected(const Select& selected) {
	std::vector<phy::Nucleus> out;
	for(const auto& [nuc, nitems, _] : selected) {
		for(u32 i=1; i <= nitems; ++i) {
			out.push_back(nuc);
		}
	}
	return out;
}
