#pragma once

/* Make sure you also build with:
 *
 * export PKG_CONFIG_PATH="$PWD/ceres-solver/build/lib/pkgconfig:$PKG_CONFIG_PATH"
 * -Iceres-solver/include -I./ -Iceres-solver/third_party/abseil-cpp -Iceres-solver/build/include -std=c++17 
 * -Lceres-solver/build/lib -Wl,rpath,ceres-solver/build/lib -lceres 
 *    $(pkg-config --libs absl_log absl_log_internal_check_op absl_strings absl_base)
 *
 *  A lot of internal dependencies. \,(>.>),/
 */

#include "Eigen/Dense"
#include "ceres/ceres.h"
#include "ceres/covariance.h"

struct GaussianCostFunctor {
	GaussianCostFunctor(double x, double y) : x_(x), y_(y) {}

	template <typename T>
	bool operator()(const T* const a0,
			const T* const a1,
			const T* const a2,
			T* residual) const {
		const T d = T(x_) - a1[0];
		const T s = a2[0];
		const T model = a0[0] * ceres::exp(-(d*d) / (T(2) * s*s));
		residual[0] = model - T(y_);
		return true;
	}

private:
	// Observations for a sample.
	const double x_;
	const double y_;
};

struct GaussianFitResult {
	int n; 
	Eigen::Vector3d params;   // [a0, a1, a2]
	Eigen::Vector3d sigma;    // 1-sigma uncertainties (std dev), from covariance diag
	Eigen::Matrix3d cov;      // full covariance (optionally scaled)
	ceres::Solver::Summary summary;

	inline double Chi2() const noexcept { return 2 * summary.final_cost; }
	inline double Chi2R() const noexcept { return Chi2() / (n - 3); };
};

GaussianFitResult FitGaussianWithUncertainty (
	const Eigen::VectorXd& x,
	const Eigen::VectorXd& y,
	Eigen::Vector3d init, // initial guess [a0,a1,a2]
	bool scale_by_residual_variance = true
) {
	if(x.size() != y.size()) throw std::runtime_error("x and y size mismatch");

	const int n = static_cast<int>(x.size());
	const int m = 3;
	if(n < m) throw std::runtime_error("Need at least 3 points for 3 parameters");
	double a0 = init[0], a1 = init[1], a2 = init[2];

	ceres::Problem problem;
	problem.AddParameterBlock(&a0, 1);
	problem.AddParameterBlock(&a1, 1);
	problem.AddParameterBlock(&a2, 1);
	
	// Enforce sigma > 0
	problem.SetParameterLowerBound(&a2, 0, 1e-12);

	for (int i = 0; i < n; ++i) {
		auto* cost = new ceres::AutoDiffCostFunction<GaussianCostFunctor, 1, 1, 1, 1> (
			new GaussianCostFunctor(x[i], y[i])
		);
		problem.AddResidualBlock(cost, nullptr, &a0, &a1, &a2);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 50;
	options.linear_solver_type = ceres::DENSE_QR; // AI overlords say this is gucci enough.
	options.minimizer_progress_to_stdout = false;

	GaussianFitResult out;
	out.n = 0;
	ceres::Solve(options, &problem, &out.summary);

	out.params << a0, a1, a2;

	// Covariance at the solution
	ceres::Covariance::Options cov_opts;
	ceres::Covariance covariance(cov_opts);

	std::vector<std::pair<const double*, const double*>> blocks = {
		{&a0, &a0}, {&a0, &a1}, {&a0, &a2},
		{&a1, &a1}, {&a1, &a2},
		{&a2, &a2}
	};

	if (!covariance.Compute(blocks, &problem)) {
		throw std::runtime_error("ceres::Covariance::Compute failed");
	}

	Eigen::Matrix3d C;
	double C00, C01, C02, C11, C12, C22;
	covariance.GetCovarianceBlock(&a0, &a0, &C00);
	covariance.GetCovarianceBlock(&a0, &a1, &C01);
	covariance.GetCovarianceBlock(&a0, &a2, &C02);
	covariance.GetCovarianceBlock(&a1, &a1, &C11);
	covariance.GetCovarianceBlock(&a1, &a2, &C12);
	covariance.GetCovarianceBlock(&a2, &a2, &C22);

	C << C00, C01, C02,
	     C01, C11, C12,
	     C02, C12, C22;

	// Optional scaling by residual variance estimate s^2 = RSS/(n-m)
	if (scale_by_residual_variance) {
		const double final_cost = out.summary.final_cost; // = 0.5 * RSS
		const int dof = n - m;
		if (dof > 0) {
			const double s2 = (2.0 * final_cost) / double(dof);
			C *= s2;
		}
	}

	out.cov = C;

	// Diagonal -> std deviations (guard against tiny negative due to floating pt errors)
	out.sigma[0] = std::sqrt(std::max(0.0, C(0,0)));
	out.sigma[1] = std::sqrt(std::max(0.0, C(1,1)));
	out.sigma[2] = std::sqrt(std::max(0.0, C(2,2)));

	return out;
}
