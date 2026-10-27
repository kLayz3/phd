#pragma once

#include "monad/accessories/primitives.hxx"
#include <limits>

namespace mnd {

static_assert(sizeof(f64) == 8, "Must hold for the implementation to be valid.");
static_assert(sizeof(f32) == 4, "Must hold for the implementation to be valid.");
static_assert(std::numeric_limits<f64>::is_iec559);
static_assert(std::numeric_limits<f32>::is_iec559);

static_assert(std::numeric_limits<f64>::radix == 2);
static_assert(std::numeric_limits<f64>::digits == 53);
static_assert(std::numeric_limits<f64>::max_exponent == 1024);

static_assert(std::numeric_limits<f32>::radix == 2);
static_assert(std::numeric_limits<f32>::digits == 24);
static_assert(std::numeric_limits<f32>::max_exponent == 128);

}
