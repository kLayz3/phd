#pragma once
#include <cstdint>
#include <cmath>

namespace phy {
namespace si {
constexpr double e = 1.602176634e-19;   // C
constexpr double c = 299792458.;        // m/s 
constexpr double me = 9.1093837139e-31; // kg
constexpr double u = 1.66053906892e-27; // kg

}

namespace nuc {
constexpr double e = 1.0;
constexpr double c = 1.0;
constexpr double u = 931.49410372;   // MeV/c^2
constexpr double me = 0.51099895069; // MeV/c^2
constexpr double Tm = 299.792458;    // MeV/(e*c) : conversion for tesla-meter
}
}


namespace phy {

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

}

