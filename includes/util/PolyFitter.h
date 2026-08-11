#pragma once

#include "../Eigen.h"
#include <algorithm>
#include <array>
#include <vector>
#include <cassert>

/* This struct serves as a fitting container.
 * for polynomials of arbitrary rank.
 * Static in a sense that the `x[N]` points remains invariant.
 * This initialization is at constructor call location.
 * `R` is the rank of the polynomial (linear => R==1 )*/

template<std::size_t N, std::size_t R>
struct StaticPolyFitter {
	static_assert(R > 0, "Fit rank (polynomial) must be greater than 0");

	using QR = Eigen::ColPivHouseholderQR<Eigen::Matrix<double, N, R+1>>;

	StaticPolyFitter() = delete;
	
	explicit StaticPolyFitter(const std::array<double, N>& x) {
		/* Eigen has bunch of expression templates, and also they vectorize very nicely.
		 * Columns are contiguous, so to preserve (kinda) cache alignment, use their API. */
		Eigen::Map<const Eigen::Matrix<double, N, 1>> xv(x.data());
		
		Eigen::Matrix<double, N, R+1> A {};
		A.col(0).setOnes();
		for(size_t i = 1; i <= R; ++i)
			A.col(i) = A.col(i-1).cwiseProduct(xv);
		qr.compute(A);
	}

	explicit StaticPolyFitter(const std::vector<double>& x) {
		assert(x.size() == N && "Vector `x` must have exactly N elements.");
		Eigen::Map<const Eigen::Matrix<double, N, 1>> xv(x.data());

		Eigen::Matrix<double, N, R+1> A {};
		A.col(0).setOnes();
		for(size_t i = 1; i <= R; ++i)
			A.col(i) = A.col(i-1).cwiseProduct(xv);
		qr.compute(A);
	}

	std::array<double, R+1> Fit(const std::array<double, N>& _y) const {
		Eigen::Map<const Eigen::Matrix<double, N, 1>> ys( _y.data() );
		
		Eigen::Matrix<double, R+1, 1> fit = qr.solve(ys);

		std::array<double, R+1> result {};
		std::copy_n(fit.data(), R+1, result.data());

		return result;
	}
	std::array<double, R+1> Fit(const std::vector<double>& _y) const {
		assert(_y.size() >= N && "Vector size must be N or bigger.");
		Eigen::Map<const Eigen::Matrix<double, N, 1>> ys( _y.data() );
		
		Eigen::Matrix<double, R+1, 1> fit =	qr.solve(ys);

		std::array<double, R+1> result;
		std::copy_n(fit.data(), R+1, result.data());

		return result;
	}

private:
	QR qr;
};

namespace mnd { namespace detail {

template<std::size_t R>
inline void PolyFit_(
	const Eigen::Ref<const Eigen::VectorXd>& x, 
	const Eigen::Ref<const Eigen::VectorXd>& y, 
	std::size_t N, std::array<double, R+1>& result
) {
	Eigen::MatrixXd A(N, R+1);
	A.col(0).setOnes();
	for(size_t i = 1; i <= R; ++i)
		A.col(i) = A.col(i-1).cwiseProduct(x);

	Eigen::Matrix<double, R+1, 1> fit = A.colPivHouseholderQr().solve(y);
	std::copy_n(fit.data(), R+1, result.data());
}

} // namespace detail
} // namespace mnd

template<std::size_t R>
void PolyFit(const std::vector<double>& x, const std::vector<double>& y, std::array<double, R+1>& result) {
	static_assert(R > 0, "Fit rank (polynomial) must be greater than 0");
	assert(((void)("Vectors must be equally sized"), x.size() == y.size()));
	
	const std::size_t N = x.size();
	assert(N >= R + 1 && "Need at least R+1 points");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), N);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), N);

	mnd::detail::PolyFit_<R>(xv, yv, N, result);
}
template<std::size_t R, std::size_t N>
void PolyFit (
	const std::array<double, N>& x, 
	const std::array<double, N>& y, 
	std::array<double, R+1>& result,
	size_t L = N
) {
	static_assert(R > 0, "Fit rank (polynomial) must be greater than 0");
	static_assert(N >= R + 1, "Static array needs to be sized least R+1.");
	assert((L >= R + 1 && L <= N) && "Need at least R+1 points, and points supplied mustn't be more than array size.");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), L);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), L);

	mnd::detail::PolyFit_<R>(xv, yv, L, result);
}

/* Fit a dataset (xi,yi) represented by two vectors `x` and `y` by a polynomial. Returns an array of coefficients. */
template<std::size_t R>
std::array<double, R+1> PolyFit(const std::vector<double>& x, const std::vector<double>& y) {
	std::array<double, R+1> res;
	PolyFit<R>(x,y, res);
	return res;
}

/* Cache friendlier version.
 * Fit a dataset (xi,yi) represented by two equally-sized arrays `x` and `y` by a polynomial. 
 * Also (optionally) supply the number of points (a slice size parameter) 
 * Returns an array of coefficients. */
template<std::size_t R, std::size_t N>
std::array<double, R+1> PolyFit(const std::array<double, N>& x, const std::array<double, N>& y, size_t L = N) {
	std::array<double, R+1> res;
	PolyFit<R>(x,y, res, L);
	return res;
}

/* Weighted least squares. 
 * If each point (xi,yi) also has a wi>0 value attached, then just by rescaling
 * A_ij' = sqrt(wi) * A_ij
 * yi' = sqrt(wi)*yi
 * we come back to ordinary least-square method. */
template<std::size_t R>
void PolyFit (
	const std::vector<double>& x, 
	const std::vector<double>& y,
	const std::vector<double>& w,
	std::array<double, R+1>& result
) {
	assert(((void)("Vectors `x` and `y` must be equally sized"), x.size() == y.size()));
	assert(((void)("Vectors `x` and `w` must be equally sized"), x.size() == w.size()));

	const std::size_t N = x.size();
	assert(N >= R + 1 && "Need at least R+1 points");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), N);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), N);
	Eigen::Map<const Eigen::VectorXd> wv(w.data(), N);

	assert((wv.array() >= 0.0).all() && "Weights must be non-negative");
	assert(wv.array().isFinite().all() && "Weights must be finite");

	const Eigen::VectorXd sqrtw = wv.array().sqrt();
	Eigen::MatrixXd A(N, R+1);
	A.col(0) = sqrtw;
	for(size_t i = 1; i <= R; ++i)
		A.col(i) = A.col(i-1).cwiseProduct(xv);

	const Eigen::VectorXd b = yv.cwiseProduct(sqrtw);

	const Eigen::Matrix<double, R+1, 1> fit = A.colPivHouseholderQr().solve(b);
	std::copy_n(fit.data(), R+1, result.data());
}

template<std::size_t R>
std::array<double, R+1> PolyFit (
	const std::vector<double>& x, 
	const std::vector<double>& y,
	const std::vector<double>& w
) {
	std::array<double, R+1> res;
	PolyFit<R>(x,y,w,res);
	return res;
}

/* Horner's algorithm: https://en.wikipedia.org/wiki/Horner%27s_method 
 * Written recursive to unroll everything. */
namespace poly {
	namespace detail {
	template<std::size_t I, std::size_t R>
	double EvalImpl__(double x, const std::array<double,R>& a) noexcept {
		if constexpr(I == R - 1)
			return a[I];
		else
			return a[I] + x * EvalImpl__<I + 1>(x, a);
	}
	inline double EvalImpl__(double x, const double* a, const int N) noexcept {
		if(N == 1) 
			return a[0];
		else
			return a[0] + x * EvalImpl__(x, a+1, N-1);
	}
	}

	template<std::size_t R>
	double Eval(const double x, const std::array<double, R>& a) noexcept {
		if constexpr(R == 0) 
			return 0.0;
		else
			return detail::EvalImpl__<0>(x, a);
	}

	inline double Eval(const double x, const std::vector<double>& a) noexcept {
		if(a.size() == 0) return 0;
		return detail::EvalImpl__(x, a.data(), static_cast<int>(a.size()) );
	}
} // namespace poly

/* Check the .cxx file to explain the reason behind explicit instantiations. */
extern template void PolyFit< 1>(const std::vector<double>& , const std::vector<double>& , std::array<double,  2>& );
extern template void PolyFit< 2>(const std::vector<double>& , const std::vector<double>& , std::array<double,  3>& );
extern template void PolyFit< 3>(const std::vector<double>& , const std::vector<double>& , std::array<double,  4>& );
extern template void PolyFit< 4>(const std::vector<double>& , const std::vector<double>& , std::array<double,  5>& );
extern template void PolyFit< 5>(const std::vector<double>& , const std::vector<double>& , std::array<double,  6>& );
extern template void PolyFit< 6>(const std::vector<double>& , const std::vector<double>& , std::array<double,  7>& );
extern template void PolyFit< 7>(const std::vector<double>& , const std::vector<double>& , std::array<double,  8>& );
extern template void PolyFit< 8>(const std::vector<double>& , const std::vector<double>& , std::array<double,  9>& );
extern template void PolyFit< 9>(const std::vector<double>& , const std::vector<double>& , std::array<double, 10>& );
extern template void PolyFit<10>(const std::vector<double>& , const std::vector<double>& , std::array<double, 11>& );
extern template void PolyFit<11>(const std::vector<double>& , const std::vector<double>& , std::array<double, 12>& );
extern template void PolyFit<12>(const std::vector<double>& , const std::vector<double>& , std::array<double, 13>& );

extern template std::array<double,  2> PolyFit< 1>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  3> PolyFit< 2>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  4> PolyFit< 3>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  5> PolyFit< 4>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  6> PolyFit< 5>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  7> PolyFit< 6>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  8> PolyFit< 7>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  9> PolyFit< 8>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 10> PolyFit< 9>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 11> PolyFit<10>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 12> PolyFit<11>(const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 13> PolyFit<12>(const std::vector<double>& , const std::vector<double>& );

extern template std::array<double,  2> PolyFit< 1>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  3> PolyFit< 2>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  4> PolyFit< 3>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  5> PolyFit< 4>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  6> PolyFit< 5>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  7> PolyFit< 6>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  8> PolyFit< 7>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double,  9> PolyFit< 8>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 10> PolyFit< 9>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 11> PolyFit<10>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 12> PolyFit<11>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
extern template std::array<double, 13> PolyFit<12>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );

extern template std::array<double, 2> PolyFit<1,2>(const std::array<double, 2>&, const std::array<double, 2>&, size_t );
extern template std::array<double, 2> PolyFit<1,3>(const std::array<double, 3>&, const std::array<double, 3>&, size_t );
extern template std::array<double, 2> PolyFit<1,4>(const std::array<double, 4>&, const std::array<double, 4>&, size_t );

/* Runtime polynomial degree version.. */
/* Input arguments to the runtime fitter should be a span type, not simply a vector.
 * But since C++17/20 issues, it's not completely trivial to do it without tagging
 * along the entire monad dependency.
 * Namely, above its matching only vector type, since the API for statically
 * known both the poly degree AND container size is the `StaticPolyFitter` */
#if !defined(MND_INCLUDE_SPAN_IS_DEFINED)
#define MND_INCLUDE_SPAN_IS_DEFINED
#if __cplusplus >= 202000L 
#	include <span>
	namespace mnd {
		template<typename T>
		using span = std::span<T>;
	}
#elif __has_include("boost/beast/core/span.hpp")
#	include "boost/beast/core/span.hpp"
	namespace mnd {
		    template<typename T>
			using span = boost::beast::span<T>;
	}
#else
#	error "Neither C++20 given, nor boost span library found. Cannot proceed."
#endif
#endif // MND_INCLUDE_SPAN_IS_DEFINED

void                PolyFit(size_t, mnd::span<const double> , mnd::span<const double> , std::vector<double>& );
std::vector<double> PolyFit(size_t, mnd::span<const double> , mnd::span<const double> );

/* Runtime polynomial degree version with weight coefficients... */
void                PolyFit(size_t, mnd::span<const double> , mnd::span<const double> , mnd::span<const double> , std::vector<double>& );
std::vector<double> PolyFit(size_t, mnd::span<const double> , mnd::span<const double> , mnd::span<const double> );

/* One extra algorithm to solve general 2D linear problem,
 * arising from solving rotational measurements:
 * x' = cos(t)*x + sin(t)*y  , AKA:
 * x' = a*x + b*y      where a^2 + b^2 == 1
 * Where `(x,y)` is the 'true' referent measurement,
 * and `x'` is what the detector gives us.
 * We want to solve for `a` and `b`.
 * We are given sequences of events: `(x,y, x')` here
 * given as the vectors `x0`, `y0`, `x` respectively.
 */

struct AngleFitResult {
	enum class Direction { X, Y };
	double a; // cos(tx)
	double b; // sin(tx)
	double Angle(const Direction ) const noexcept; // tx
};

struct AngleOffsetFitResult {
	AngleFitResult t;
	double c; // -dx*cos(t.tx) - dy*sin(t.tx)
};

AngleFitResult FitAngle (
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
);

AngleOffsetFitResult FitAngleOffset (
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
);
