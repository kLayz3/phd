#pragma once

/* Eigen's fancy AVX512 SIMD instructors give a stroke
 * to GCC 14.2.0, so hide the standard diagnostics so that the tty
 * can live in peace. */

#if defined(__GNUC__) && !defined(__clang__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#	pragma GCC diagnostic ignored "-Wunused-variable"
#elif defined(__clang__)
#	pragma clang diagnostic push
#	pragma clang diagnostic ignored "-Wunused-variable"
#	pragma clang diagnostic ignored "-Wunused-variable"
#endif

#include "Eigen/Dense"

#if defined(__GNUC__) && !defined(__clang__)
#	pragma GCC diagnostic pop
#elif defined(__clang__)
#	pragma clang diagnostic pop
#endif
