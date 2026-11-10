#pragma once

#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

#include <vector>
#include <type_traits>

#include "monad/accessories/typedefs.hxx"

namespace mnd {

template<typename T>
inline static constexpr T nan = std::numeric_limits<T>::quiet_NaN();

template<typename T>
using Interpolator = boost::math::interpolators::pchip<std::vector<T>>;

namespace detail {

template<typename T>
struct Bounds {
	T lo_x, hi_x;
	T lo_y, hi_y;

};

} // namespace detail

namespace bounds {
enum MatchDerivative : u32 {
	Left  = 1 << 0,
	Right = 1 << 1,
};
} // namespace bounds

/* Generally we interpolate inbetween via Cubic Hermite, and extrapolate outside
 * by some fits or ad-hoc functions. */
template <
	typename T,
	typename LoExtrapolator,
	typename HiExtrapolator
> struct Fitter : private detail::Bounds<T>, private Interpolator<T> {
	static_assert(std::is_invocable_r_v<T, const LoExtrapolator&, T>,
		"Lower edge extrapolator must be invocable as: auto f(T) -> T");
	static_assert(std::is_invocable_r_v<T, const HiExtrapolator&, T>,
		"High edge extrapolator must be invocable as: auto f(T) -> T");
	
	using Ip = Interpolator<T>;
	using Bn = detail::Bounds<T>;

	Fitter() = delete;

	/* In general, the extrapolators can be lambdas (closure objects) which can hold
	 * complex internal state, such as fit parameters. Now, if taking the extrapolators
	 * by reference, would be quite hard to verify the lifetimes (Rust fixes this btw bruh).
	 * As such, we take the extrapolators *by-value*, which is negligible for function-pointers 
	 * or small-ish lambdas. */

	/* The ctor aborts on empty vector. Hehe. In full release mode, first access if undefined, but does not segfault
	 * (unless vec capacity is exactly 0). But then the ctor afterwards panics. So all in all - we are in the clear.
	 *
	 * The value of h(hi_x) will be made to match the interpolator's value at the edge, which is y.back().
	 * I.e.: h(x.back())  == y.back()
	 *       l(x.front()) == y.front()
	 *  */
	explicit Fitter(std::vector<T> x, std::vector<T> y, LoExtrapolator l, HiExtrapolator h,
		u32 match_derivative = bounds::Left | bounds::Right) :
		Bn(checked_bounds(x,y)),
		Ip(std::move(x), std::move(y),
			get_derivative(match_derivative & bounds::Left , l, Bn::lo_x),
			get_derivative(match_derivative & bounds::Right, h, Bn::hi_x)
		),
		ext_lo(std::move(l)), ext_hi(std::move(h)),
		ext_lo_offset( Bn::lo_y - ext_lo(Bn::lo_x) ),
		ext_hi_offset( Bn::hi_y - ext_hi(Bn::hi_x) )
		{}

	T operator()(T value) const noexcept {
		// Firstly, check if value is a NAN. Give it instantly back then.
		// Infinities get handled by their extrapolators
		if(std::isnan(value))
			return value;

		if(value < Bn::lo_x) {
			return ext_lo(value) + ext_lo_offset;
		}
		else if(value > Bn::hi_x) {
			return ext_hi(value) + ext_hi_offset;
		}

		return Ip::operator()(value);
	}

	template<typename F>
	static T get_derivative(bool match, F&& functor, T x) {
		if(!match) // Let PCHIP choose its default slope.
			return nan<T>;
		return boost::math::differentiation::
			finite_difference_derivative(functor, x);
	}

	static Bn checked_bounds(const std::vector<T>& x, const std::vector<T>& y) {
		assert(x.size() >= 4);
		assert(x.size() == y.size());

		return {x.front(), x.back(), y.front(), y.back()};
	}

private:
	LoExtrapolator ext_lo;
	HiExtrapolator ext_hi;

	T ext_lo_offset;
	T ext_hi_offset;
};

template <
	typename T,
	typename LoExtrapolator,
	typename HiExtrapolator
> Fitter(::std::vector<T>, ::std::vector<T>, LoExtrapolator, HiExtrapolator, u32 = 0)
	-> Fitter<T,LoExtrapolator,HiExtrapolator>;

template<typename T>
struct PowerLaw {
	T amplitude;
	T exponent;

	inline T operator()(T E) const noexcept {
		return amplitude * std::pow(E, exponent);
	}

	PowerLaw inverse() const noexcept {
		const double one_over_exponent = 1.0/exponent;
		return {
			.amplitude = std::pow(this->amplitude, -one_over_exponent),
			.exponent = one_over_exponent
		};
	}
};

template<typename T>
PowerLaw(T, T) -> PowerLaw<T>;


template<typename T>
struct InvalidFitter {
	T operator()(T ) const noexcept { return nan<T>; }
};

} // namespace mnd
