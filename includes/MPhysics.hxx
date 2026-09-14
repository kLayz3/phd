#pragma once

#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <string_view>
#include <ostream>

#include "util/CLI.h"
#include "util/Option.hxx"
#include "util/FromChars.h"
#include "StaticNuclearTable.hh"

namespace phy {
namespace si {

constexpr double e = 1.602176634e-19;   // C
constexpr double c = 299792458.;        // m/s
constexpr double me = 9.1093837139e-31; // kg
constexpr double u = 1.66053906892e-27; // kg

} // namespace si

namespace nuc {

constexpr double e = 1.0;
constexpr double c = 1.0;
constexpr double u = 931.49410372;   // MeV/c^2
constexpr double me = 0.51099895069; // MeV/c^2
constexpr double be = 13.6e-6;       // MeV, binding energy of electron in 1H atom
constexpr double Tm = 299.792458;    // MeV/(e*c) , conversion for tesla-meter

} // namespace nuc

namespace detail {
/* Allow doing runtime asserts for calulcations.
 * Most API otherwise just returns a quiet NAN on failure. */
inline constexpr bool debug_ = false;
} // namespace detail

template<typename T>
using Option = mnd::Option<T>;

struct Brho_t {
	double value;
};
struct Beta_t {
	double value;
};
struct EKin_t {
	double value; // Kinetic energy per nucleon.
};

inline constexpr char chem_symbols[][3] = {
	"",
	"H",  "He",
	"Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
	"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
	"K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
	"Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
	"Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
	"Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
	"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
	"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
	"Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
	"Tl", "Pb", "Bi", "Po", "At", "Rn",
	"Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
	"Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
	"Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
	"Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
constexpr uint32_t n_chem_symbols = (uint32_t)sizeof(chem_symbols) / sizeof(*chem_symbols);

inline constexpr mnd::Option<uint32_t> GetZ(std::string_view symbol) noexcept {
	const char first  = symbol[0];
	const char second = symbol.size() >= 2 ? symbol[1] : '\0';

	/* Even tho strlen(chem_symbol[i]) can be 1, it's safe to access chem_symbol[1]
	 * since it will be just the null-terminator. */
	for(uint32_t z = 1; z < n_chem_symbols; ++z)
		if(chem_symbols[z][0] == first && chem_symbols[z][1] == second)
			return mnd::Some{z};

	return mnd::None;
}

/* Return back atomic number based on the `string_view` object provided.
 * Will mutate the view and move it one past, based on how many (1 or 2)
 * characters it munched to deduce the chemical element. */
inline constexpr mnd::Option<uint32_t> Z(std::string_view& symbol) noexcept {
	if(symbol.size() >= 2) {
		/* Try to match symbols with 2 characters. */
		if(const auto z = GetZ(symbol.substr(0,2)); z.is_some()) {
			symbol.remove_prefix(2);
			return z;
		}
	}

	if(!symbol.empty()) {
		/* Try to match symbols with just 1 char. */
		if(const auto z = GetZ(symbol.substr(0,1)); z.is_some()) {
			symbol.remove_prefix(1);
			return z;
		}
	}

	return mnd::None;
}

struct Nucleus {
	uint32_t A = 0;
	uint32_t Z = 0;

	/* Number of electrons, possibly not given. */
	mnd::Option<uint16_t> N_electrons;

	/* Returns the mass of the nucleus. */
	inline double mass() const noexcept {
		const double m0 = ::phy::mass(A,Z);
		if constexpr(detail::debug_) {
			if(!std::isfinite(m0))
				MND_THROW("Mass value looked up as invalid value. A=%u, Z=%u", A,Z);
		}
		return m0 + N_electrons.map([this](const auto value) {
				return value * nuc::me -
					(this->Z) * (this->Z) * nuc::be *
					( std::min((int)value, 2) + std::max((int)value - 2, 0) / 4.0 );
			}).value_or(0.0);
	}
	inline double mass_per_nucleon() const noexcept {
		return mass() / A;
	}

	inline Option<uint32_t> charge_state() const noexcept {
		return N_electrons.and_then([this](uint16_t ne) -> Option<uint32_t> {
			if(ne > Z) return mnd::None;
			return mnd::Some{Z - ne};
		});
	}

	friend inline std::ostream& operator<<(std::ostream& os, Nucleus const& n) {
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
};

namespace cvt {

inline constexpr double get_beta_from_betagamma2(double beta_gamma2) noexcept {
	return sqrt( beta_gamma2 / (1 + beta_gamma2) );
}
inline constexpr double get_beta_from_betagamma(double beta_gamma) noexcept {
	return beta_gamma / sqrt(1 + beta_gamma * beta_gamma);
}

}; // namespace cvt

inline double Gamma(double b) noexcept {
	if constexpr(detail::debug_) {
		if(b >= 1 || b < 0)
			MND_THROW("Requested beta velocity %.2f out of bounds.", b);
	}
	if(b >= 1) return std::numeric_limits<double>::infinity();
	if(b < 0)  return NAN;
	return 1.0 / std::sqrt(1.0 - b*b);
}

/* Brho is a number corresponding to units of tesla-meter. */
inline double BetaGamma(uint32_t A, uint32_t Q, Brho_t brho) noexcept {
	const double m = mass(A,Q);
	return brho.value * Q * nuc::Tm / m;
}
inline double BetaGamma2(uint32_t A, uint32_t Q, EKin_t e) noexcept {
	const double m = mass(A,Q);
	const double tmp = (e.value * A)/ m + 1;
	return tmp*tmp - 1;
}
/* Sometimes A,Q is known at comptime. */
template<uint32_t A, uint32_t Q>
constexpr double BetaGamma(Brho_t brho) noexcept {
	constexpr double m = mass<A,Q>();
	constexpr double q_tm_over_m = Q * nuc::Tm / m;
	return brho.value * q_tm_over_m;
}
template<uint32_t A, uint32_t Q>
constexpr double BetaGamma2(EKin_t e) noexcept {
	constexpr double m = mass<A,Q>();
	const double tmp = (e.value * A)/ m + 1;
	return tmp*tmp - 1;
}

/* Brho is a number corresponding to units of tesla-meter. */
inline double Beta(uint32_t A, uint32_t Q, Brho_t brho) noexcept {
	const double beta_gamma = BetaGamma(A,Q,brho);
	return cvt::get_beta_from_betagamma(beta_gamma);
}
inline double Beta(uint32_t A, uint32_t Q, EKin_t e) noexcept {
	const double beta_gamma2 = BetaGamma2(A,Q,e);
	return cvt::get_beta_from_betagamma2(beta_gamma2);
}

/* Sometimes A,Q is known at comptime. */
template<uint32_t A, uint32_t Q>
constexpr double Beta(Brho_t brho) noexcept {
	const double beta_gamma = BetaGamma<A,Q>(brho);
	return cvt::get_beta_from_betagamma(beta_gamma);
}
template<uint32_t A, uint32_t Q>
constexpr double Beta(EKin_t e) noexcept {
	const double beta_gamma2 = BetaGamma2<A,Q>(e);
	return cvt::get_beta_from_betagamma2(beta_gamma2);
}

/* Kinetic energy per nucleon, if (A,Q,beta) are known. */
inline double EKin(uint32_t A, uint32_t Q, Beta_t beta) noexcept {
	const double g = Gamma(beta.value);
	const double gb = g * beta.value;
	const double frac = std::max(sqrt(gb*gb + 1.0) - 1.0, 0.0);
	return frac * mass(A,Q) / A;
}
/* Kinetic energy per nucleon, if (A,Q,brho) are known. */
inline double EKin(uint32_t A, uint32_t Q, Brho_t brho) noexcept {
	const double gb = BetaGamma(A, Q, brho);
	const double frac = std::max(sqrt(gb*gb + 1.0) - 1.0, 0.0);
	return frac * mass(A,Q) / A;
}

template<uint32_t A, uint32_t Q>
constexpr double EKin(Beta_t beta) noexcept {
	constexpr double mass_per_nucleon = mass<A,Q> / A;
	const double g = Gamma(beta.value);
	const double gb = g * beta.value;
	const double frac = std::max(sqrt(gb*gb + 1.0) - 1.0, 0.0);
	return frac * mass_per_nucleon;
}
template<uint32_t A, uint32_t Q>
constexpr double EKin(Brho_t brho) noexcept {
	constexpr double mass_per_nucleon = mass<A,Q> / A;
	const double gb = BetaGamma<A,Q>(brho.value);
	const double frac = std::max(sqrt(gb*gb + 1.0) - 1.0, 0.0);
	return frac * mass_per_nucleon;
}
/* ^^^ Here we make an overload with strongly-typed beta vs. brho and not give
 * a generic `double` overload. This was a source of headache. \_(>_<)_/ */

inline double EKin(const Nucleus& n, Beta_t beta) {
	const double g = Gamma(beta.value);
	const double gb = g * beta.value;
	const double frac = std::max(sqrt(gb*gb + 1.0) - 1.0, 0.0);
	return frac * n.mass_per_nucleon();
}

namespace literals {

inline constexpr Brho_t operator""_brho(long double value) noexcept {
	return { static_cast<double>(value) };
}

} // namespace literals
} // namespace phy
