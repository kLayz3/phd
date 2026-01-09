#pragma once

#include "Eigen/Dense"
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


template<size_t R>
void PolyFit(const std::vector<double>& x, const std::vector<double>& y, std::array<double, R+1>& result) {
	static_assert(R >= 0, "Fit rank (polynomial) must be greater than 0");

	assert(((void)("Vectors must be equally sized"), x.size() == y.size()));
	
	const std::size_t N = x.size();
	assert(N >= R + 1 && "Need at least R+1 points");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), N);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), N);

	Eigen::MatrixXd A(N, R+1);
	A.col(0).setOnes();
	for(size_t i = 1; i <= R; ++i)
		A.col(i) = A.col(i-1).cwiseProduct(xv);

	Eigen::Matrix<double, R+1, 1> fit = A.colPivHouseholderQr().solve(yv);

	std::copy_n(fit.data(), R+1, result.data());
}

template<size_t R>
std::array<double, R+1> PolyFit(const std::vector<double>& x, const std::vector<double>& y) {
	std::array<double, R+1> res;
	PolyFit<R>(x,y, res);
	return res;
}
