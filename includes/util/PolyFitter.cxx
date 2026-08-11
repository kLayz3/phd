#include "PolyFitter.h"

/* This hack is used to instantiate the template for small-ish N's.
 * This has the advantage that then these calls can then be used in ROOT macros.
 *
 * Problem is that Eigen is heavily optimized from rank 2 onward for algorithms such
 * as colPivHouseholderQr() decomposition and ROOT's cling when it tries to compile the source of macro,
 * it sees opaque Eigen optimizations (SIMD stuff), and you will catch the weirdest segfault in your life.
 * Who would've thought that making a C++ interpreter would make users' lives easier?
 *
 * Don't ask me why, but I was swearing in >4 different languages debugging this... */

template void PolyFit< 1>(const std::vector<double>& , const std::vector<double>& , std::array<double,  2>& );
template void PolyFit< 2>(const std::vector<double>& , const std::vector<double>& , std::array<double,  3>& );
template void PolyFit< 3>(const std::vector<double>& , const std::vector<double>& , std::array<double,  4>& );
template void PolyFit< 4>(const std::vector<double>& , const std::vector<double>& , std::array<double,  5>& );
template void PolyFit< 5>(const std::vector<double>& , const std::vector<double>& , std::array<double,  6>& );
template void PolyFit< 6>(const std::vector<double>& , const std::vector<double>& , std::array<double,  7>& );
template void PolyFit< 7>(const std::vector<double>& , const std::vector<double>& , std::array<double,  8>& );
template void PolyFit< 8>(const std::vector<double>& , const std::vector<double>& , std::array<double,  9>& );
template void PolyFit< 9>(const std::vector<double>& , const std::vector<double>& , std::array<double, 10>& );
template void PolyFit<10>(const std::vector<double>& , const std::vector<double>& , std::array<double, 11>& );
template void PolyFit<11>(const std::vector<double>& , const std::vector<double>& , std::array<double, 12>& );
template void PolyFit<12>(const std::vector<double>& , const std::vector<double>& , std::array<double, 13>& );

template std::array<double,  2> PolyFit< 1>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  3> PolyFit< 2>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  4> PolyFit< 3>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  5> PolyFit< 4>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  6> PolyFit< 5>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  7> PolyFit< 6>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  8> PolyFit< 7>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  9> PolyFit< 8>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 10> PolyFit< 9>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 11> PolyFit<10>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 12> PolyFit<11>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 13> PolyFit<12>(const std::vector<double>& , const std::vector<double>& );

template std::array<double,  2> PolyFit< 1>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  3> PolyFit< 2>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  4> PolyFit< 3>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  5> PolyFit< 4>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  6> PolyFit< 5>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  7> PolyFit< 6>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  8> PolyFit< 7>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double,  9> PolyFit< 8>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double, 10> PolyFit< 9>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double, 11> PolyFit<10>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double, 12> PolyFit<11>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );
template std::array<double, 13> PolyFit<12>(const std::vector<double>& , const std::vector<double>& , const std::vector<double>& );

template std::array<double, 2> PolyFit<1,2>(const std::array<double, 2>&, const std::array<double, 2>&, size_t );
template std::array<double, 2> PolyFit<1,3>(const std::array<double, 3>&, const std::array<double, 3>&, size_t );
template std::array<double, 2> PolyFit<1,4>(const std::array<double, 4>&, const std::array<double, 4>&, size_t );

void PolyFit (
	size_t R,
	mnd::span<const double> x,
	mnd::span<const double> y,
	std::vector<double>& result
) {
	assert(R >= 1 && "Runtime polynomial fit: rank passed must be >= 1."); 
	assert(((void)("Vectors must be equally sized"), x.size() == y.size()));
	
	const size_t Npts = x.size();
	assert(Npts >= R + 1 && "Runtime polynomial fit: need at least R+1 points.");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), Npts);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), Npts);
	Eigen::MatrixXd A(Npts, R+1);
	A.col(0).setOnes();
	for(size_t i = 1; i <= R; ++i)
		A.col(i) = A.col(i-1).cwiseProduct(xv);

	const Eigen::VectorXd fit = A.colPivHouseholderQr().solve(yv); 
	result.assign(fit.data(), fit.data() + fit.size());
}
[[ nodiscard ]]
std::vector<double> PolyFit (
	size_t N,
	mnd::span<const double> x, 
	mnd::span<const double> y
) {
	std::vector<double> res;
	PolyFit(N, x, y, res);
	return res;
}

void PolyFit (
	size_t R, 
	mnd::span<const double> x, 
	mnd::span<const double> y,
	mnd::span<const double> w,
	std::vector<double>& result
) {
	assert(R >= 1 && "Runtime polynomial fit: rank passed must be >= 1."); 
	assert(((void)("Vectors must be equally sized"), x.size() == y.size()));
	assert(((void)("Vectors `x` and `w` must be equally sized"), x.size() == w.size()));
	
	const size_t Npts = x.size();
	assert(Npts >= R + 1 && "Runtime polynomial fit: need at least R+1 points.");

	Eigen::Map<const Eigen::VectorXd> xv(x.data(), Npts);
	Eigen::Map<const Eigen::VectorXd> yv(y.data(), Npts);
	Eigen::Map<const Eigen::VectorXd> wv(w.data(), Npts);
	
	assert((wv.array() >= 0.0).all() && "Weights must be non-negative");
	assert(wv.array().isFinite().all() && "Weights must be finite");
	const Eigen::VectorXd sqrtw = wv.array().sqrt();

	Eigen::MatrixXd A(Npts, R+1);
	A.col(0) = sqrtw;
	for(size_t i = 1; i <= R; ++i)
		A.col(i) = A.col(i-1).cwiseProduct(xv);
	
	const Eigen::VectorXd b = yv.cwiseProduct(sqrtw);

	Eigen::VectorXd fit = A.colPivHouseholderQr().solve(b); 
	result.assign(fit.data(), fit.data() + fit.size());
}
[[ nodiscard ]]
std::vector<double> PolyFit (
	size_t R,
	mnd::span<const double> x,
	mnd::span<const double> y,
	mnd::span<const double> w
) {
	std::vector<double> res;
	PolyFit(R, x, y, w, res);
	return res;
}

double AngleFitResult::Angle(AngleFitResult::Direction d) const noexcept {
	switch(d) {
		case Direction::X: 
			return std::atan2(b, a);
		case Direction::Y:
			return -M_PI / 2 + std::atan2(b, a);
	}
	return 0;
}

AngleFitResult FitAngle (
    const std::vector<double>& x0, // Referent measurement (along the device's axis)
    const std::vector<double>& y0, // Referent measurement (orthogonal to the device's axis)
    const std::vector<double>& x   // Device's measurements
) {
	const std::size_t N = x.size();
	assert(x0.size() == N && "FitAngle(): vectors `x0` and `x` must be identically sized.");
	assert(y0.size() == N && "FitAngle(): vectors `y0` and `x` must be identically sized.");
	assert(N >= 2 && "Must supply more than 2 points.");

	Eigen::MatrixXd A(N, 2);
	Eigen::VectorXd b(N);
	for (std::size_t i = 0; i < N; ++i) {
		A(i, 0) = x0[i];
		A(i, 1) = y0[i];
		b(i)    = x[i];
	}
	Eigen::Vector2d fit = A.colPivHouseholderQr().solve(b);

	return {
		fit(0), // ~ cos(tx)
		fit(1)  // ~ sin(tx)
	};
}

AngleOffsetFitResult FitAngleOffset (
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
) {
	const std::size_t N = x.size();
	assert(x0.size() == N && "FitAngleOffset(): vectors `x0` must be identically sized.");
	assert(y0.size() == N && "FitAngleOffset(): vectors `y0` must be identically sized.");
	assert(N >= 3 && "Must supply more than 2 points.");

	Eigen::MatrixXd A(N, 3);
	Eigen::VectorXd rhs(N);
	for (std::size_t i = 0; i < N; ++i) {
		A(i, 0) = x0[i];
		A(i, 1) = y0[i];
		A(i, 2) = 1.0;
		rhs(i)  = x[i];
	}
	Eigen::Vector3d fit = A.colPivHouseholderQr().solve(rhs);

	AngleOffsetFitResult out;
	out.t.a = fit(0);
	out.t.b = fit(1);
	out.c   = fit(2);
	return out;
}; 
