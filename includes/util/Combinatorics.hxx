#pragma once

#include <cstdint>
#include <array>
#include <cassert>

namespace mnd::combi {

inline constexpr
std::size_t n_choose_k(uint32_t n, uint32_t k) {
    if(k > n) return 0;
    if(k * 2 > n) k = n-k;
    if(k == 0) return 1;

    std::size_t result = static_cast<std::size_t>(n);
    for(uint32_t i = 2; i <= k; ++i ) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

/* Extract the i-th unique combination (subset) of size M.
 * For a range of size N, there are N_C_M unique combinations whose indices 
 * {a(0), a(1), ..., a(M-1)} can be ordered lexicographically, and then
 * represented as a bitmask.
 * E.g for n=5, k=3, there are n_choose_k combinations == 10
 *
 *  [0] -> 00111  {0,1,2}
 *  [1] -> 01011  {0,1,3}
 *  [2] -> 10011  {0,1,4}
 *  [3] -> 01101  {0,2,3}
 *  [4] -> 10101  {0,2,4}
 *  [5] -> 11001  {0,3,4}
 *  [6] -> 01110  {1,2,3}
 *  [7] -> 10110  {1,2,4}
 *  [8] -> 11010  {1,3,4}
 *  [9] -> 11100  {2,3,4}
 */
inline constexpr uint32_t combination (
	uint32_t n,
	uint32_t k,
	std::size_t index
) {
	assert(n <= 32 && "mnd::combi::combination(), sequence size must be <= 32");
	assert(k <= n && "mnd::combi::combination(), subset size > sequence size ?");
    assert(index < n_choose_k(n, k) && "mnd::combi::combination(), combi index > n combinations ?");

	/* FYI this is a nice Leetcode problem. Check LC #77 :-) */
	uint32_t mask = 0;
    uint32_t first = 0;

    for(uint32_t pos = 0; pos < k; ++pos) {
        const uint32_t remaining = k - pos - 1;
        for(uint32_t i = first; i < n; ++i) {
            const std::size_t count = n_choose_k(n - i - 1, remaining);
            if(index < count) {
                mask |= uint32_t{1} << i;
                first = i + 1;
                break;
            }
            index -= count;
        }
    }
	return mask;
}

template<uint32_t N, uint32_t K>
constexpr auto make_combo_lookup_table() {
	static_assert(N <= 32);
	static_assert(K <= N);

	std::array<uint32_t, n_choose_k(N, K)> result {};

	for(std::size_t i = 0; i < result.size(); ++i)
		result[i] = combination(N, K, i);

	return result;
}

template<uint32_t N, uint32_t K>
inline constexpr auto combo_lookup_table =
	make_combo_lookup_table<N, K>();

} // namespace mnd::combi
