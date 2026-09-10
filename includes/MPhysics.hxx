#pragma once

#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <string_view>

#include "util/Option.hxx"

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
constexpr double be = 13.6e-6;       // MeV
constexpr double Tm = 299.792458;    // MeV/(e*c) : conversion for tesla-meter

} // namespace nuc

inline double Gamma(double b) noexcept {
	if(b >= 1) return std::numeric_limits<double>::infinity();
	if(b < 0)  return NAN;
	return 1.0 / std::sqrt(1.0 - b*b);
}

/* Brho is a number corresponding to units of tesla-meter. */
inline double Beta(uint32_t Q, uint32_t A, double brho) noexcept {
	/* Alpha = Q/A * (e * Tm / (uc) ) * brho */
	constexpr double e_times_tm_over_uc = 1.0 * phy::nuc::Tm / (phy::nuc::u * 1.0);
	double alpha = Q*brho / A * e_times_tm_over_uc;
	return alpha / sqrt(1 + alpha*alpha); 
}
/* Sometimes Q is known. */
template<uint32_t Q>
double Beta(uint32_t A, double brho) noexcept {
	constexpr double Q_times_e_times_tm_over_uc = Q * 1.0 * phy::nuc::Tm / (phy::nuc::u * 1.0);
	double alpha = brho / A * Q_times_e_times_tm_over_uc;
	return alpha / sqrt(1 + alpha*alpha); 
}
/* Sometimes Q,A is known. */
template<uint32_t Q, uint32_t A>
double Beta(double brho) noexcept {
	constexpr double Q_over_A_times_e_times_tm_over_uc = static_cast<double>(Q)/A * 1.0 * phy::nuc::Tm / (phy::nuc::u * 1.0);
	double alpha = brho * Q_over_A_times_e_times_tm_over_uc;
	return alpha / sqrt(1 + alpha*alpha); 
}

/* Kinetic energy per nucleon. */
inline double EKin(double beta) noexcept {
	double g = phy::Gamma(beta);
	double gb = g*beta;
	double frac = sqrt(gb*gb + 1.0) - 1;
	if(frac < 0) frac = 0.;
	return frac * nuc::u; 
}
/* Kinetic energy per nucleon, if (Q,A,brho) are known. */
inline double EKin(uint32_t Q, uint32_t A, double brho) noexcept {
	double beta = phy::Beta(Q,A,brho);
	return EKin(beta);
}
/* Sometimes Q is known. */
template<uint32_t Q>
double EKin(uint32_t A, double brho) noexcept {
	double beta = phy::Beta<Q>(A,brho);
	return EKin(beta);
}
/* Sometimes Q,A is known. */
template<uint32_t Q, uint32_t A>
double EKin(double brho) noexcept {
	double beta = phy::Beta<Q,A>(brho);
	return EKin(beta);
}

constexpr mnd::Option<uint32_t> GetZ(std::string_view symbol) noexcept {
	constexpr char symbols[][3] = {
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

	const char first  = symbol[0];
	const char second = symbol.size() >= 2 ? symbol[1] : '\0';

	for(uint32_t z = 1; z < sizeof(symbols)/ sizeof(*symbols); ++z)
		if(symbols[z][0] == first && symbols[z][1] == second)
			return mnd::Some{z};

	return mnd::None;
}

/* Return back atomic number based on the `string_view` object provided.
 * Will mutate the view and move it one past, based on how many (1 or 2) 
 * characters it munched to deduce the chemical element. */
constexpr mnd::Option<uint32_t> Z(std::string_view& symbol) noexcept {
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

} //namespace phy

#include "StaticNuclearTable.hh"

namespace phy {

struct Nucleus {
	uint32_t Z, A;

	/* Number of electrons */
	mnd::Option<uint16_t> N_electrons;

	inline double mass() const noexcept {
		return ::phy::mass(Z,A)
			+ N_electrons.map([this](const auto value) {
				return value * nuc::me +
					(this->Z) * (this->Z) * nuc::be *
					( std::min((int)value, 2) + std::max((int)value - 2, 0) / 4.0 );
			}).value_or(0.0);
	}
};

} //namespace phy
