#pragma once

#include <array>

/* In the experiment, there are a few things we can assume and "leave" to the comptime.
 * E.g., we cannot directly measure the velocity/kin.energy of the fragment directly *before* the S2 target.
 * But from sims, LISE++, can roughly describe it to have minor corrections to different measurements. */

namespace mnd::assume {
namespace s2 {

/* Kinetic energy loss [AMeV] from the entrance of S2 to the Be target. */
inline constexpr double loss_upto_target = 11.3;

/* Kinetic energy loss [AMeV] from the entrance of S2 Be target up to the exit. */
inline constexpr double loss_in_target = 59.0;

/* 9Be target density. */
inline constexpr double target_density = 1.845; // [g/cm^3]

/* S2 is dispersive, so the incoming energy before an element slightly depends on the x- position.
 * dE = f(x), where f(0) = 0. This linear interpolation is eye-balled from LISE++ , minding the
 * fact we use Go4 convention of left-right, which is reversed compared to LISE! */
inline constexpr auto e_dispersion = std::array{
	0.0,          // [MeV] ; 0 by definition
	+14.0 / 80.0  // [MeV / mm]
};

/* Sometime vertex calculation could blow up. In this case,
 * just clamp the value into this `x_bound`. */
inline constexpr auto x_bound = std::array{
	-40.0, // [mm]
	40.0   // [mm]
};

} // namespace s2
} // namespace mnd::assume
