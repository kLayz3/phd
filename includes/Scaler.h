#pragma once

#include "monad/monad.hxx"

template<uint32_t N>
struct Scaler {
	static_assert(N > 1 && N <= 32, "Template parameter for `Scaler` must be >1 && <= 32.");
	static constexpr u32 mask_ = static_cast<u32>((1ull << N) - 1);
	static constexpr i64 wrap_point_ = 1ll << (N-1);

	u64 cumulative {};

	u32 prev_data {};
	u32 curr_data {};
	u32 increment {};

	bool initialized_ {false};

	Scaler() = default;
	
	inline void assign(u32 fresh) noexcept {
		if constexpr(N != 32)
			fresh &= mask_;
		if(!initialized_) {
			curr_data = fresh;
			initialized_ = true;
			return;
		}
		prev_data = curr_data;
		curr_data = fresh;
		increment = calc_increment();
		cumulative += static_cast<uint64_t>(increment);
	}

	u32 calc_increment() const noexcept {
		if(curr_data >= prev_data) {
			return curr_data - prev_data;
		}
		/* Possible miscounting! */
		if(prev_data - curr_data < static_cast<u32>(wrap_point_)) {
			WARN("Backwards counting in scaler struct. Prev = %u, curr = %u\n", prev_data, curr_data);
			return 0;
		}
		/* Wrap-around. */
		return (u32)((1ll << N) + (i64)curr_data - (i64)prev_data);
	}

	/* `x` and `y` should be at most one wrap-around different.
	 * Wraparound is recognized if their distance is more than half of the scale. */
	static i32 calc_diff(u32 x, u32 y) noexcept {
		x &= mask_; y &= mask_;

		i64 raw_diff = static_cast<i64>(x) - static_cast<i64>(y);
		if(raw_diff > wrap_point_) { // `y` is one wrap ahead.
			return static_cast<i32>(raw_diff - (1ll<<N));
		}
		else if(raw_diff < -wrap_point_) { // `x` is one wrap ahead.
			return static_cast<i32>(raw_diff + (1ll<<N));
		}
		else {
			return static_cast<i32>(raw_diff);
		}
	}

	virtual ~Scaler() = default;
	ClassDef(Scaler, 1);
};

