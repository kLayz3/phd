#pragma once

#include "Eigen/Dense"
#include <algorithm>
#include <array>

/* This struct serves as a fitting container.
 * for polynomials of arbitrary rank. */

template<std::size_t N, size_t R>
struct StaticPolyFitter {
	static_assert(R > 0, "Fit rank (polynomial) must be greater than 0");
	StaticPolyFitter() = default;
	StaticPolyFitter(std::array<double, N> const& x) {
		for(size_t i = 0; i < N; ++i) {
			double powx = 1;
			A(i, 0) = 1.0;
			for(size_t r = 1; r <= R; ++r) {
				powx *= x[i];
				A(i,r) = powx; 
			}
		}
	}
	StaticPolyFitter(std::vector<double> const& x) {
		assert(x.size() == N && "Vector unequally sized?");

		for(size_t i = 0; i < N; ++i) {
			double powx = 1;
			A(i, 0) = 1.0;
			for(size_t r = 1; r <= R; ++r) {
				powx *= x[i];
				A(i,r) = powx; 
			}
		}
	}

	std::array<double, R+1> Fit(const std::array<double, N>& _y) {
		Eigen::Map<const Eigen::Matrix<double, N, 1>> ys( _y.data() );
		
		Eigen::Matrix<double, R+1, 1> fit = A.colPivHouseholderQr().solve(ys);

		std::array<double, R+1> result;
		std::copy_n(fit.data(), R+1, result);

		return result;
	}
	std::array<double, R+1> Fit(const std::vector<double>& _y) {
		assert(_y.size() >= N && "Vector size must be N or bigger.");
		Eigen::Map<const Eigen::Matrix<double, N, 1>> ys( _y.data() );
		
		Eigen::Matrix<double, R+1, 1> fit = A.colPivHouseholderQr().solve(ys);

		std::array<double, R+1> result;
		std::copy_n(fit.data(), R+1, result.data());

		return result;
	}

private:
	Eigen::Matrix<double, N, R+1> A;
};
