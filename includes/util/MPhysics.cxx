#include "util/MPhysics.h"
#include "util/Geometry.h"
#include "monad/monad.hxx"

#include <algorithm>
#include <cassert>
#include <numeric>

using namespace phy;

bool Nucleus::valid() const noexcept {
	return A > 0 &&
	       Z > 0 &&
	       std::isfinite(mass()) &&
	      (N_electrons.is_none() ||
	       N_electrons.unwrap() <= Z);
}

double Nucleus::mass() const noexcept {
	const double m0 = ::phy::mass(A,Z);
	if constexpr(detail::debug_) {
		if(!std::isfinite(m0))
			MND_THROW("Mass value looked up as invalid value. A=%u, Z=%u", A,Z);
	}
	return m0 +
		N_electrons.map([this](const auto value) {
			return value * nuc::me -
				(this->Z) * (this->Z) * nuc::be *
				( std::min((int)value, 2) + std::max((int)value - 2, 0) / 4.0 );
		}).value_or(0.0);
}
double Nucleus::mass_per_nucleon() const noexcept {
	return mass() / A;
}

Option<uint32_t> Nucleus::charge_state() const noexcept {
	return N_electrons.and_then([this](uint16_t ne) -> Option<uint32_t> {
		if(ne > Z) return mnd::None;
		return mnd::Some{Z - ne};
	});
}

Option<AtomicNumber> Nucleus::ChemElem() const noexcept {
	return ToAtomicNumber(Z);
}
std::string Nucleus::to_string() const noexcept {
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}

std::string Nucleus::chem_to_string(
	Represent repr,
	bool add_nucleon_number
) const noexcept {
	return ChemElem()
		.map([this, repr, add_nucleon_number](AtomicNumber z) -> std::string {
			if(z == AtomicNumber::He && A == 4) {
				if(repr == Represent::Rootex)
					return "#alpha";
				else if(repr == Represent::Pythex)
					return R"(\alpha)";
				else
					return "α";
			}
			if(z == AtomicNumber::H && A == 1) { // proton
				if(repr == Represent::Pythex)
					return R"(\mathrm{p})";
				else
					return "d";
			}
			if(z == AtomicNumber::H && A == 2) { // deuteron
				if(repr == Represent::Pythex)
					return R"(\mathrm{d})";
				else
					return "d";
			}

			std::string_view symbol = magic_enum::enum_name(z);
			std::string extra = (add_nucleon_number
				? mnd::sstrcat("{}^{", mnd::utos(A), '}')
				: "");

			if(repr == Represent::Rootex) // "{}^{7}Be" | "Be"
				return std::move(extra) + std::string{symbol};
			else if(repr == Represent::Pythex) // "\mathrm{{}^{7}Be}" | "\mathrm{Be}"
				return mnd::sstrcat(
					R"(\mathrm{)", std::move(extra), symbol, '}');
			else
				return std::string{symbol};
		})
		.value_or("");
}

Nucleus Nucleus::get_ion(std::string_view str,
	bool bad_parse_is_error,
	const char* desc
) {
	if(!desc) ERROR("Description string is a nullptr? Pass \"\" empty string instead!\n");
	std::string_view view = mnd::trim(str);

	Option<u32> A = mnd::stou_munch(view);
	if(A.is_none()) {
		auto msg = mnd::msg("mnd::RunsheetState::get_ion(desc: \"%s\") for input string: \'%*s\', "
			"couldn't parse out the mass number (A). Got `None`. "
			"Setting the value to default: %u\n", desc, (int)view.length(), view.data(), DEFAULT_A_PRIMARY);
		if(!bad_parse_is_error) {
			WARN("%s", msg);
		} else {
			throw std::invalid_argument(msg);
		}
		A = mnd::Some{DEFAULT_A_PRIMARY};
	}

	Option<u32> Z = phy::Z(view);
	if(Z.is_none()) {
		auto msg = mnd::msg("mnd::RunsheetState::get_ion(desc: \"%*s\") for input string: \'%s\' "
			"(current view: \'%*s\'), "
			"couldn't parse out the atomic number (Z). Got `None`. "
			"Setting the value to default: %u\n", desc, (int)str.length(), str.data(),
			(int)view.length(), view.data(), DEFAULT_A_PRIMARY);
		if(!bad_parse_is_error) {
			WARN("%s", msg);
		} else {
			throw std::invalid_argument(msg);
		}
		Z = mnd::Some{DEFAULT_Z_PRIMARY};
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
			return (Z.unwrap() >= q) ? mnd::Some{(u16)(Z.unwrap() - q)} : mnd::Some{u16{0}};
		})
	};
	/* There can be more junk after the parse, but that's fine. */
}

bool phy::operator==(Nucleus const& lhs, Nucleus const& rhs) noexcept {
	return lhs.A == rhs.A &&
	       lhs.Z == rhs.Z &&
	       lhs.N_electrons == rhs.N_electrons;
}
bool phy::operator!=(Nucleus const& lhs, Nucleus const& rhs) noexcept {
	return !(lhs == rhs);
}

std::ostream& phy::operator<<(std::ostream& os, Nucleus const& n) {
	if(n.A == 0 || n.Z == 0)
		return os << "{}";
	if(n.Z >= n_chem_symbols)
		return os << "{invalid Z}";
	
	// E.g. 56Fe26+
	return os << n.A << chem_symbols[n.Z]
		<< n.charge_state().map([](uint32_t q) -> std::string {
			return mnd::utos(q) + "+";
		}).value_or("");
}

std::istream& phy::operator>>(std::istream& is, Nucleus& nucleus) {
	std::string token;
	if(is >> token)
		nucleus = Nucleus::get_ion(token, true, "phy::operator>>");

	return is;
}

Nucleus phy::operator&(Nucleus const& lhs, Nucleus const& rhs) {
	Option<u16> ne = lhs.N_electrons.
		and_then([&rhs](u16 value) -> Option<u16> {
			if(rhs.N_electrons.is_some())
				return mnd::Some{(u16)(value + rhs.N_electrons.unwrap())};
				/* Thanks C++ for having a retarded-ass rule that
				* uint16_t + uint16_t gives int. ROFLMAO.
				* If I'm too stupid to care for wrap-arounds then let me be stupid,
				* not have this idiotic rule of IMPLICIT promotions on binary ops.
				* Please God, get me to Rust I can't this anymore. */
			else
				return mnd::None;
		});

	return Nucleus {
		.A = lhs.A + rhs.A,
		.Z = lhs.Z + rhs.Z,
		.N_electrons = ne
	};
}

Nucleus& Nucleus::operator&=(Nucleus const& rhs) {
	*this = *this & rhs;
	return *this;
}

template<>
double phy::rho<phy::RhoExpressionType::full>(
	mnd::span<const Nucleus> daughters,
	mnd::span<const mnd::geom::Line3D> tracks
) {
	if constexpr(detail::debug_) {
		assert(daughters.size() == tracks.size() && daughters.size() > 1 &&
			"phy::rho(): nuclei span and tracks are unequal lengths?");
		[[maybe_unused]] bool valid_masses = std::all_of(daughters.cbegin(), daughters.cend(),
			[](const Nucleus& n) { return n.valid(); });
		assert(valid_masses &&
			"phy::rho(): nuclei span seems to have invalid daughter (can't take its mass)?");
	}
	
	const u32 N = (u32)daughters.size();
	const double M = std::accumulate(daughters.begin(), daughters.end(), 0.0,
		[](double sum, const auto& n) { return sum + n.mass(); });

	double rho2 = 0; /* [mrad]^2 */
	for(u32 i=0; i<N; ++i) {
		for(u32 j=i+1; j<N; ++j) {
			double mfactor =
				daughters[i].mass() * daughters[j].mass() / (M * mp);

			double theta = MRAD_CVT * tracks[i].AngleRelativeTo(tracks[j]);
			
			rho2 += mfactor * theta * theta;
		}
	}

	return std::sqrt(rho2);
}

template<>
double phy::rho<phy::RhoExpressionType::p>(
	mnd::span<const Nucleus> daughters,
	mnd::span<const mnd::geom::Line3D> tracks
) {
	if constexpr(detail::debug_) {
		assert(daughters.size() == tracks.size() && daughters.size() > 1 &&
			"phy::rho(): nuclei span and tracks are unequal lengths or length 0?");
		[[maybe_unused]] bool valid_masses = std::all_of(daughters.cbegin(), daughters.cend(),
			[](const Nucleus& n) { return n.valid(); });
		assert(valid_masses &&
			"phy::rho(): nuclei span seems to have invalid daughter (can't take its mass)?");
	}

	/* As per RNFOOTTrack convention. Heavy Ion (HI) must be at the front.
	 * The validity of the views we don't check here. Must be verified by the caller. */
	const u32 N = (u32)daughters.size() - 1; /* Number of protons. */
	const double Mh = daughters.front().mass();
	const double mass_frac = mp / Mh; /* < 1.0 */

	const mnd::geom::Line3D& heavy_track = tracks.front();

	double rho2 = 0; /* [mrad]^2 */
	for(u32 i=1; i <= N; ++i) {
		double theta_i = MRAD_CVT * heavy_track.AngleRelativeTo(tracks[i]);
		rho2 += theta_i * theta_i;
	}
	
	{
		double rho2_pp = 0;
		for(u32 i=1; i <= N; ++i) {
			for(u32 j=i+1; j <= N; ++j) {
				double theta_ij = MRAD_CVT * tracks[i].AngleRelativeTo(tracks[j]);
				rho2_pp += theta_ij * theta_ij;
			}
		}
		rho2_pp *= mass_frac;
		rho2 += rho2_pp;
	} // leave rho2_pp out of scope, to not mess up.

	return std::sqrt( rho2 / (1 + N * mass_frac) );
}

template<>
double phy::rho<phy::RhoExpressionType::equinuclear>(
	mnd::span<const Nucleus> daughters,
	mnd::span<const mnd::geom::Line3D> tracks
) {
	if constexpr(detail::debug_) {
		assert(daughters.size() == tracks.size() && daughters.size() > 1 &&
			"phy::rho(): nuclei span and tracks are unequal lengths or length 0?");
		[[maybe_unused]] bool valid_masses = std::all_of(daughters.cbegin(), daughters.cend(),
			[](const Nucleus& n) { return n.valid(); });
		assert(valid_masses &&
			"phy::rho(): nuclei span seems to have invalid daughter (can't take its mass)?");
	}
	
	/* The validity of the views we don't check here. Must be verified by the caller. */
	const u32 N = (u32)daughters.size(); /* Number of daughters. */
	const double M0 = daughters.back().mass();
	
	double rho2 = 0;
	for(u32 i=0; i < N; ++i) {
		for(u32 j=i+1; j < N; ++j) {
			double theta_ij = MRAD_CVT * tracks[i].AngleRelativeTo(tracks[j]);
			rho2 += theta_ij * theta_ij;
		}
	}

	return std::sqrt(
		M0 / (N*mp) * rho2
	);
}

template<>
double phy::rho<phy::RhoExpressionType::unknown>(
	mnd::span<const Nucleus> daughters,
	mnd::span<const mnd::geom::Line3D> tracks
) {
	(void)daughters, (void)tracks;
	MND_THROW("Bad enum argument given to the template. Try again.");
}

double phy::rho(
	RhoExpressionType e,
	mnd::span<const Nucleus> daughters,
	mnd::span<const mnd::geom::Line3D> tracks
) {
	switch(e) {
		case(RhoExpressionType::p):
			return rho<RhoExpressionType::p>(daughters, tracks);
		case(RhoExpressionType::equinuclear):
			return rho<RhoExpressionType::equinuclear>(daughters, tracks);
		case(RhoExpressionType::full):
			return rho<RhoExpressionType::full>(daughters, tracks);
		case(RhoExpressionType::unknown):
			return rho<RhoExpressionType::unknown>(daughters, tracks);
		default:
			mnd::unreachable();
	}
}

