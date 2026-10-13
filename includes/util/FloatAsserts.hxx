#pragma once
#include <cstdint>
#include <limits>

namespace mnd {

using f64 = double;
using f32 = float;

using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 =  int16_t;
using i32 =  int32_t;
using i64 =  int64_t;

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
