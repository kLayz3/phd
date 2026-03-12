#pragma once

#include "../Eigen/Dense"
#include <algorithm>
#include <array>
#include <vector>
#include <cassert>

/* This struct serves as a fitting container.
 * for polynomials of arbitrary rank.
 * Static in a sense that the `x[N]` points remains invariant. 
 * `R` is the rank of the polynomial (linear => R==1 )*/

template<std::size_t N, std::size_t R>
struct StaticPolyFitter {
	static_assert(R >= 0, "Fit rank (polynomial) must be greater than 0");

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

template<std::size_t R>
void PolyFit(const std::vector<double>& x, const std::vector<double>& y, std::array<double, R+1>& result) {
	static_assert(R >= 0, "Fit rank (polynomial) must be greater than 0");

	assert(((void)("Vectors must be equally sized"), x.size() == y.size()));
	
	const std::size_t N = x.size();
	assert(N >= R + 1 && "Need at least R+1 points");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), N);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), N);

	PolyFit_<R>(xv, yv, N, result);
}

template<std::size_t R>
std::array<double, R+1> PolyFit(const std::vector<double>& x, const std::vector<double>& y) {
	std::array<double, R+1> res;
	PolyFit<R>(x,y, res);
	return res;
}

/* Weighted least squares. 
 * If each point (xi,yi) also has a wi>0 value attached, then just by rescaling
 * xi' = sqrt(wi)*xi
 * yi' = sqrt(wi)*yi
 * we come back to ordinary least-square method. */

/* Weights vector `w` must have all elements >= 0. 
 * This check **is not** reinforced at runtime. */
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

	Eigen::VectorXd x_reduced = wv.array().sqrt() * xv.array();
	Eigen::VectorXd y_reduced = wv.array().sqrt() * yv.array();
	
	PolyFit_<R>(x_reduced, y_reduced, N, result);
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

	template<std::size_t R>
	double Eval(const double x, const std::array<double, R>& a) noexcept {
		if constexpr(R == 0) 
			return 0.0;
		else
			return EvalImpl__<0>(x, a);
	}

	inline double Eval(const double x, const std::vector<double>& a) noexcept {
		if(a.size() == 0) return 0;
		return EvalImpl__(x, a.data(), static_cast<int>(a.size()) );
	}
}
