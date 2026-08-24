#pragma once

/* MND defined way to manipulate quiet NAN's and embed information
 * opaquely into the mantissa. */

/* Cool info: I use `memcpy` to represent the bytes of a double, and then copy it back.
 * Reason is that the hack: `u64 value = *(u64*)(&val);` is infact UB under C++'s strict aliasing rule:
 * https://gist.github.com/shafik/848ae25ee209f698763cffee272a58f8
 * Namely, compiler can assume that u64* pointer cannot mutate the original double, and can completely optimize away
 * (e.g., ignore the store) any write to that adress through the aliased `u64*` pointer.
 * In general, this is only really UB if binding the repr to a function return, e.g.:
 * __attribute__((noinline)) u64* as_bits(double* p) {
 *    return reinterpret_cast<std::uint64_t*>(p);
 * }
 * double f() { f64 d = 1.0; u64_t* repr = as_bits(&d); *repr = 0; return d; }
 * This gives different `f()` value if compiled under `-fno-strict-aliasing` or without!
 * No clue why GCC even under `-Wextra` doesn't warn about this. */

#include <cstdint>
#include <cstring>
#include <limits>

namespace mnd {

using f64 = double;
using f32 = float;

using u32 = uint32_t;
using u64 = uint64_t;

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

/* Construct a quiet NAN with a specific mantissa payload
 * given as bytes of an 8-byte value. */
inline f64 make_nan64(u64 payload) noexcept {
	const u64 repr = (0x7ff8000000000000ULL
		| (payload & 0x0007ffffffffffffULL));

	f64 result;
	std::memcpy(&result, &repr, sizeof repr);
	return result;
}
inline f32 make_nan32(u32 payload) noexcept {
	const u32 repr = (0x7fc00000U
		| (payload & 0x003fffffU));

	f32 result;
	std::memcpy(&result, &repr, sizeof repr);
	return result;
}

inline u64 get_nan_payload(f64 value) noexcept {
	u64 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & 0x0007ffffffffffffULL;
}
inline u32 get_nan_payload(f32 value) noexcept {
	u32 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & 0x003fffffU;
}
template<typename T>
auto get_nan_payload(T ) = delete;

/* Extract bits from a floating point value indicated by the bitmask argument. */
inline u64 get_fbits(f64 value, u64 bitmask) noexcept {
	u64 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & bitmask;
}
inline u32 get_fbits(f32 value, u32 bitmask) noexcept {
	u32 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & bitmask;
}
template<typename T, typename U>
auto get_fbits(T , U ) = delete;

template<u32 index>
bool get_bit(f64 value) noexcept {
	static_assert(index < 64);
	constexpr u64 mask = (u64{1} << index);
	u64 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & mask;
}
template<u32 index>
bool get_bit(f32 value) noexcept {
	static_assert(index < 32);
	constexpr u32 mask = (u32{1} << index);
	u32 repr;
	std::memcpy(&repr, &value, sizeof repr);
	return repr & mask;
}
template<u32 , typename T>
auto get_bit(T ) = delete;

/* Index is the bit index, counted from [31:0] for f32
 * and [63:0] for f64 */
template<u32 index, bool val = true>
void set_bit(f64& value) noexcept {
	static_assert(index < 64);
	constexpr u64 mask = (u64{1} << index);
	u64 repr;
	std::memcpy(&repr, &value, sizeof repr);
	if constexpr(val) {
		repr |= mask;
	} else {
		repr &= ~mask;
	}
	std::memcpy(&value, &repr, sizeof repr);
}
template<u32 index, bool val = true>
void set_bit(f32& value) noexcept {
	static_assert(index < 32);
	constexpr u32 mask = (u32{1} << index);
	u32 repr;
	std::memcpy(&repr, &value, sizeof repr);
	if constexpr(val) {
		repr |= mask;
	} else {
		repr &= ~mask;
	}
	std::memcpy(&value, &repr, sizeof repr);
}

/* Set the lowest bits: [nbits-1:0] to have the properly
 * truncated binary representation of an integer value `value`. */
template<u32 nbits>
void set_low_bits(f64& number, u64 value) noexcept {
	static_assert(nbits > 0 && nbits <= 52, "Must request to set between 1 and 52 bits.");
	u64 repr;
	std::memcpy(&repr, &number, sizeof repr);
	constexpr u64 BIT_MASK  = (u64{1} << nbits) - 1;
	constexpr u64 BIT_NMASK = ~BIT_MASK;
	
	repr = (repr & BIT_NMASK) | (value & BIT_MASK);
	std::memcpy(&number, &repr, sizeof repr);
}
template<u32 nbits>
void set_low_bits(f32& number, u32 value) noexcept {
	static_assert(nbits > 0 && nbits <= 23, "Must request to set between 1 and 23 bits.");
	u32 repr;
	std::memcpy(&repr, &number, sizeof repr);
	constexpr u32 BIT_MASK  = (u32{1} << nbits) - 1;
	constexpr u32 BIT_NMASK = ~BIT_MASK;
	
	repr = (repr & BIT_NMASK) | (value & BIT_MASK);
	std::memcpy(&number, &repr, sizeof repr);
}

} // namespace mnd
