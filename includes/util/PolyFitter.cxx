#include "PolyFitter.h"
#include "util/PolyFitter.hxx"

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
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& x
) {
	const std::size_t N = x.size();
	assert(x0.size() == N && "FitAngle(): vectors `x0` must be identically sized.");
	assert(y0.size() == N && "FitAngle(): vectors `y0` must be identically sized.");
	assert(N >= 2 && "Must supply more than 2 points.");

	Eigen::MatrixXd A(N, 2);
	Eigen::VectorXd b(N);
	for (std::size_t i = 0; i < N; ++i) {
		A(i, 0) = x0[i];
		A(i, 1) = y0[i];
		b(i)    = x[i];
	}
	Eigen::Vector2d fit = A.colPivHouseholderQr().solve(b);

	[[ maybe_unused ]] const double c = fit(0);  // ~ cos(tx)
	[[ maybe_unused ]] const double d = fit(1);  // ~ sin(tx)
	
	return {
		fit(0),
		fit(1)
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
