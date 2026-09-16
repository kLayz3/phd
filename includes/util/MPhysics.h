#pragma once

#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <string_view>
#include <ostream>

#include "magic_enum/magic_enum.hpp"
#include "util/Option.hxx"
#include "util/FromChars.h"
#include "util/Geometry.h"
#include "StaticNuclearTable.hh"

namespace phy {

template<typename T>
using Option = ::mnd::Option<T>;

using namespace std::string_literals;

constexpr static auto None = ::mnd::None;

namespace si {

inline constexpr double e = 1.602176634e-19;   // C
inline constexpr double c = 299792458.;        // m/s
inline constexpr double me = 9.1093837139e-31; // kg
inline constexpr double u = 1.66053906892e-27; // kg

} // namespace si

namespace nuc {

inline constexpr double e = 1.0;
inline constexpr double c = 1.0;
inline constexpr double u = 931.49410372;   // MeV/c^2
inline constexpr double me = 0.51099895069; // MeV/c^2
inline constexpr double be = 13.6e-6;       // MeV, binding energy of electron in 1H atom
inline constexpr double Tm = 299.792458;    // MeV/(e*c) , conversion for tesla-meter
inline constexpr double mp = mass<1,1>();   // MeV/c^2, proton mass

} // namespace nuc

inline constexpr double mp = nuc::mp;

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

enum class AtomicNumber : uint32_t {
	H = 1, He,
	Li, Be, B, C, N, O, F, Ne,
	Na, Mg, Al, Si, P, S, Cl, Ar,
	K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni,
	Cu, Zn, Ga, Ge, As, Se, Br, Kr,
	Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd,
	Ag, Cd, In, Sn, Sb, Te, I, Xe,
	Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd,
	Tb, Dy, Ho, Er, Tm, Yb, Lu,
	Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
	Tl, Pb, Bi, Po, At, Rn,
	Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm,
	Bk, Cf, Es, Fm, Md, No, Lr,
	Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn,
	Nh, Fl, Mc, Lv, Ts, Og
};
inline constexpr uint32_t n_chem_symbols =
	static_cast<uint32_t>(AtomicNumber::Og) + 1;

namespace detail {
struct chem_symbol_table {
	char data[n_chem_symbols][3] = {'\0'};
};

inline constexpr auto chem_symbol_storage_ = []() constexpr {
	chem_symbol_table result{};

	for(const auto& entry : magic_enum::enum_entries<AtomicNumber>()) {
		const auto z = static_cast<uint32_t>(entry.first);
		const auto name = entry.second;

		if(z == 0 || z >= n_chem_symbols || name.empty() || name.size() > 2)
			throw std::invalid_argument("Invalid chem element?");

		for(size_t i = 0; i < name.size(); ++i)
			result.data[z][i] = name[i];
	}

	return result;
}();
} // namespace detail

inline constexpr auto& chem_symbols = detail::chem_symbol_storage_.data;

inline constexpr mnd::Option<AtomicNumber> ToAtomicNumber(uint32_t value) noexcept {
	if(value < static_cast<uint32_t>(AtomicNumber::H) ||
	   value > static_cast<uint32_t>(AtomicNumber::Og))
		return mnd::None;

	return mnd::Some{static_cast<AtomicNumber>(value)};
}

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
	static constexpr uint32_t DEFAULT_A_PRIMARY = 0;
	static constexpr uint32_t DEFAULT_Z_PRIMARY = 0;

	uint32_t A = DEFAULT_A_PRIMARY;
	uint32_t Z = DEFAULT_Z_PRIMARY;

	/* Number of electrons, possibly not given. */
	Option<uint16_t> N_electrons;

	bool valid() const noexcept;

	/* Returns the mass of the nucleus. */
	double mass() const noexcept;
	double mass_per_nucleon() const noexcept;
	Option<uint32_t> charge_state() const noexcept;
	Option<AtomicNumber> ChemElem() const noexcept;
	std::string to_string() const noexcept;

	Nucleus& operator&=(Nucleus const& );
	
	enum struct Represent { Rootex = 1, Pythex = 2, Normal = 3 };

	/* Represent the nuclear charge as a chemical element, subjective to either
	 * [1]: Standard TLatex (ROOT)-like math formatting
	 * [2]: Standard Latex (Python)-like math formatting
	 * [3]: Usual format as used in text. */
	std::string chem_to_string(
		Represent repr = Represent::Normal,
		bool add_nucleon_number = false
	) const noexcept;

	/* Parse a string into the Nucleus state.
	 * get_ion("16O6+") => Nucleus{ .A = 16, .Z = 8, .N_electrons = Some{2} }
	 * get_ion("")      => Nucleus{ .A = DEFAULT_A_PRIMARY, .Z = DEFAULT_Z_PRIMARY, .N_electrons = None }
	 * , second argument converts a bad parse into an std::invalid_argument throw.
	 * Third argument is only a description string for a possible error message. */
	[[ nodiscard ]]
	static Nucleus get_ion(std::string_view str,
		bool bad_parse_is_error = false,
		const char* desc = ""
	);
};
bool operator==(Nucleus const& , Nucleus const& ) noexcept;
bool operator!=(Nucleus const& , Nucleus const& ) noexcept;
std::ostream& operator<<(std::ostream& , Nucleus const& );
std::istream& operator>>(std::istream& , Nucleus& );

/* Perform a "fusion" of two nuclei. */
Nucleus operator&(Nucleus const& , Nucleus const& );

namespace literals {
inline Nucleus operator""_n(const char* text, std::size_t size) {
	return Nucleus::get_ion(
		std::string_view{text, size},
		true,
		"operator\"\"_n"
	);
}
} // namespace literals

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

/* Calculate the ρ-value given a sequence of tracks and masses.
 * The base equation used is the general one, as used by k14y23:
 * [0]: ρ = sqrt( sum_{i=1}^N sum_{j=i+1}^N  m(i)*m(j)/(M * mp) * θ(i,j)^2 )
 * The equation used can also be specified, for either [1] N proton-decays,
 * e.g.: 8C -> 4He + p + p + p + p
 * , or [2] breakup into N identical particles,
 * e.g.: 9C -> 3He + 3He + 3He
 * [1]: ρ = 1/sqrt(1 + N*mp/mH) sqrt{ sum_{i=1}^N θ(H,pi)^2 + (mp/mH) sum_{i=1}^N sum_{j=i+1}^N θ(pi,pj)^2 }
 * [2]: ρ = sqrt{m0 / (N*mp)} * sqrt{ sum_{i=1}^N sum_{j=i+1}^N θ(i,j)^2 } */
enum class RhoExpressionType : u32 {
	full        = 0, // Reaction: X => y_1 + y_2 + ... + y_N
	p           = 1, // Reaction: X => Y + Np
	equinuclear = 2, // Reaction: X => N y
	unknown     = 3, // Ditto
};

/* Conversion factor from standard mnd::geom API (which gives [rad]) */
inline constexpr double MRAD_CVT = 1'000.0;

template<enum RhoExpressionType = RhoExpressionType::full>
double rho(mnd::span<const Nucleus> , mnd::span<const mnd::geom::Line3D> );

template<> double rho<RhoExpressionType::full>(
	mnd::span<const Nucleus> ,
	mnd::span<const mnd::geom::Line3D>
);
template<> double rho<RhoExpressionType::p>(
	mnd::span<const Nucleus> ,
	mnd::span<const mnd::geom::Line3D>
);
template<> double rho<RhoExpressionType::equinuclear>(
	mnd::span<const Nucleus> ,
	mnd::span<const mnd::geom::Line3D>
);
template<> double rho<RhoExpressionType::unknown>(
	mnd::span<const Nucleus> ,
	mnd::span<const mnd::geom::Line3D>
);
double rho(
	RhoExpressionType ,
	mnd::span<const Nucleus> ,
	mnd::span<const mnd::geom::Line3D>
); // Runtime variant.

namespace literals {

inline constexpr Brho_t operator""_brho(long double value) noexcept {
	return { static_cast<double>(value) };
}

} // namespace literals
} // namespace phy
