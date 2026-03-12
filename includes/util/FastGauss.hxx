/* This small library is a routine for fast gauss fitting.
 * When the number of (x,y) pairs going into the fit is 3/4/5 and known at compile time. 
 * It's used to make quick FOOT cluster fit around the maximum highlighted strip. 
 * We drop the ceres dependency as it's too clunky. */

/* Algorithm used is Levenberg–Marquardt algorithm (LM)
 * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm */

/* In a Gauss Fitting procedure, the optimal amplitude can be calculated analytically. As such
 * the *full* Jacobian folds down to 2 parameters not three. Therefore, conceptually the fit is split like this:
 * 
 * model(x; mu,sigma) = a0(mu,sigma) * phi(x; mu,sigma)  
 * Where: phi(mu, sigma) = exp(-(x-mu)^2 / (2*sigma^2)) */

/* In all the algorithms, small `s` refers to log(sigma), as such the following is *always* true:
 * s == log(sigma) <=> sigma = exp(s) */

#pragma once
//#define FGF_LIB_ADD_PRINTER_
//#define FGF_LIB_DEBUG_

#include "Eigen/Dense"
#include <cmath>
#include <algorithm>

#ifdef FGF_LIB_ADD_PRINTER_
#	include <iostream>
#	include <iomanip>
#endif

struct FGResult {
	// Fit parameters
	double a0 = 0.0;     // Amplitude
	double mu = 0.0;     // Mean
	double sigma = 1.0;  // Width

	// Uncertainties (sqrt of covariance diagonal)
	double a0_std = 0.0;
	double mu_std = 0.0;
	double sigma_std = 0.0;

	// Diagnostic stuff
	bool ok = false;
	int iters = 0;
	double final_cost = 0.0;  // 0.5 * sum(residual^2) after weighting
};

#ifdef FGF_LIB_ADD_PRINTER_
inline std::ostream& operator<<(std::ostream& os, const FGResult& r) {
	os << std::fixed << std::setprecision(3)
	   << "A: "   << r.a0 << "±" << r.a0_std
	   << ", μ: " << r.mu << "±" << r.mu_std 
	   << ", σ: " << r.sigma << "±" << r.sigma_std
	   << " === OK: " << r.ok << ", iters: " << r.iters
	   << ", cost: " << r.final_cost;
	return os;
}
#endif

namespace fgf {

/* Small container to park the a0 dependency in the picture. */
struct AmpAndDerivs {
  double a0 = 0.0;
  double da0_dmu = 0.0;
  double da0_ds  = 0.0;
};

/* Given mu,s compute phi, its derivatives, and analytic amplitude a0 plus its derivatives.
 * phi_i = exp(-(x-mu)^2 / (2*sigma^2) ) , which is 'normal' amplitude gaussian. Only two parameter dependent. 
 * Also cache the intermediate model calculations in the *phi* buffers. */
template<int N>
AmpAndDerivs amplitude_and_derivs (
    const Eigen::Matrix<double, N, 1>& x,
    const Eigen::Matrix<double, N, 1>& y,
    double mu, double s,
    Eigen::Matrix<double, N, 1>& phi,
    Eigen::Matrix<double, N, 1>& dphi_dmu,
    Eigen::Matrix<double, N, 1>& dphi_ds
) {
	const double sigma = std::exp(s);
	const double inv_sig2 = 1.0 / (sigma * sigma);

	for(int i = 0; i < N; ++i) {
		const double d = x[i] - mu;
		const double e = std::exp(-(d*d) * 0.5 * inv_sig2);
		phi[i] = e;

		/* dphi/dmu = phi * (d / sigma^2) */
		dphi_dmu[i] = e * (d * inv_sig2);

		/* dphi/ds  = phi * (d^2 / sigma^2) */
		dphi_ds[i]  = e * ((d*d) * inv_sig2);
	}

	/* a0 = (phiᵀ y) / (phiᵀ phi) == u/v  , vector notation. */
	const double u = phi.dot(y);
	const double v = phi.squaredNorm();

	AmpAndDerivs out;

	// guard against degenerate v
	if(v <= 0.0 || !std::isfinite(v)) return out;

	out.a0 = u / v;

	/* Derivatives:
	 * du/dp = sum (dphi/dp) y_i
	 * dv/dp = sum 2*phi_i (dphi/dp) */
	const double du_dmu = dphi_dmu.dot(y);
	const double dv_dmu = 2.0 * phi.dot(dphi_dmu);

	const double du_ds  = dphi_ds.dot(y);
	const double dv_ds  = 2.0 * phi.dot(dphi_ds);

	const double inv_v2 = 1.0 / (v * v);

	out.da0_dmu = (v * du_dmu - u * dv_dmu) * inv_v2;
	out.da0_ds  = (v * du_ds  - u * dv_ds ) * inv_v2;

	return out;
}

template<std::size_t N>
double cost_half_sqnorm(const Eigen::Matrix<double, N, 1>& r) noexcept {
	return 0.5 * r.squaredNorm();
}

template<std::size_t N>
void eval_all (
	const Eigen::Matrix<double, N, 1>& x,
	const Eigen::Matrix<double, N, 1>& y,
	double mu, double s,
	double inv_meas,
	Eigen::Matrix<double, N, 1>& phi,
	Eigen::Matrix<double, N, 1>& dphi_dmu,
	Eigen::Matrix<double, N, 1>& dphi_ds,
	AmpAndDerivs& amp,
	Eigen::Matrix<double, N, 1>& r,
	Eigen::Matrix<double, N, 2>& J
) {
	amp = fgf::amplitude_and_derivs<N>(x, y, mu, s, phi, dphi_dmu, dphi_ds);
	for(std::size_t i = 0; i < N; ++i) {
		const double f = amp.a0 * phi[i];
		r[i] = (y[i] - f) * inv_meas;
		const double dr_dmu = -(amp.da0_dmu * phi[i] + amp.a0 * dphi_dmu[i]) * inv_meas;
		const double dr_ds  = -(amp.da0_ds  * phi[i] + amp.a0 * dphi_ds[i])  * inv_meas;
		J(i, 0) = dr_dmu;
		J(i, 1) = dr_ds;
	}
}

} // namespace fgf

static constexpr int    FG_DEFAULT_MAX_ITERS =  6;
static constexpr double FG_DEFAULT_Y_SIGMA   =  0.0; 
static constexpr double FG_DEFAULT_LAMBDA0   =  1e-3;
static constexpr double FG_DEFAULT_TOL_STEP  =  1e-10;
static constexpr double FG_DEFAULT_TOL_GRAD  =  1e-12;

/* init = {mu, sigma}, must be relatively good. */
template<size_t N>
FGResult FastGaussFit (
	const Eigen::Matrix<double, N, 1>& x,
    const Eigen::Matrix<double, N, 1>& y,
    const Eigen::Vector2d& init, // {mu, sigma}
    int max_iters   = FG_DEFAULT_MAX_ITERS,
    double y_sigma  = FG_DEFAULT_Y_SIGMA,  // Optionally known noise. for FOOT: around ~2.8 ADC units.
    double lambda0  = FG_DEFAULT_LAMBDA0,  // Lambda parameter in LM algorithm.
    double tol_step = FG_DEFAULT_TOL_STEP,
    double tol_grad = FG_DEFAULT_TOL_GRAD 
) noexcept {
	using VecN  = Eigen::Matrix<double, N, 1>;
	using MatN2 = Eigen::Matrix<double, N, 2>;
	using Vec2  = Eigen::Matrix<double, 2, 1>;
	using Mat2  = Eigen::Matrix<double, 2, 2>;

	FGResult out{};

	double mu = init[0];
	double sigma0 = std::max(1e-15, init[1]);
	double s = std::log(sigma0);
	fgf::AmpAndDerivs amp, amp_t;
	
	double lambda = lambda0; 

	VecN phi, dphi_dmu, dphi_ds;
	VecN r;
	MatN2 J;

	const double inv_meas = (y_sigma > 0.0) ? (1.0 / y_sigma) : 1.0;

	/* Placeholder values. */
	fgf::eval_all<N>(x, y, mu, s, inv_meas, phi, dphi_dmu, dphi_ds, amp, r, J);
	double c = fgf::cost_half_sqnorm<N>(r);

	for(int it=0; it < max_iters; ++it) {
		out.iters = it + 1;

		const Mat2 H = J.transpose() * J;
		const Vec2 g = J.transpose() * r;
	
		if(g.norm() < tol_grad) break;
		
		// Levenberg–Marquardt step: (H + λI) δ = -g
		Mat2 A = H + lambda * Mat2::Identity();
		Vec2 delta = A.ldlt().solve(-g);

		if(!delta.allFinite()) {
			out.ok = false;
			return out;
		}
		if(delta.norm() < tol_step) break;

		/* Trial. See if new values (mu', s') = (mu,s) + delta
		 * are better than previously. */
		const double mu_t = mu + delta[0];
		const double s_t  = s  + delta[1];
		VecN r_t;
		MatN2 J_t;
		
		fgf::eval_all<N>(x, y, mu_t, s_t, inv_meas, phi, dphi_dmu, dphi_ds, amp_t, r_t, J_t);
		double c_t = fgf::cost_half_sqnorm<N>(r_t);

#ifdef FGF_LIB_DEBUG_
		std::cout << (c_t < c ? "  ACCEPT " : "  REJECT ")
          << " trial_cost=" << c_t << "\n";
#endif
		if(std::isfinite(c_t) and c_t <= c) {
			/* Accept this iteration. */
			amp = amp_t; mu = mu_t; s = s_t;
			r = r_t; J = J_t; c = c_t;
			lambda = std::max(1e-18, lambda * 0.3); // Close. Go roughly newton-like
		} else {
			lambda = std::min(1e18, lambda * 3.0); // Far. Go back to gradient descent
		}

#ifdef FGF_LIB_DEBUG_
		std::cout << "it=" << it
          << " mu=" << mu
          << " sigma=" << std::exp(s)
          << " a0=" << amp.a0
          << " cost=" << c
          << " lambda=" << lambda
          << " |g|=" << g.norm()
          << " delta={" << delta[0] << ", " << delta[1] << "}\n";
#endif
	}
	
	out.a0 = amp.a0;
	out.mu = mu;
	out.sigma = std::exp(s);
	out.final_cost = c;

	// Covariance: Cov([mu,s]) ≈ scale * (JᵀJ)^-1
	Mat2 H = J.transpose() * J;
	Eigen::LDLT<Mat2> ldlt(H);
	if(ldlt.info() != Eigen::Success) {
		out.ok = false;
		return out;
	}
	Mat2 Hinv = ldlt.solve(Mat2::Identity());

	const int dof = N - 2;
	double scale = 1.0;
	if(y_sigma > 0.0) {
		scale = 1.0; // already normalized residuals => absolute scale
	} else if(dof > 0) {
		const double rss = 2.0 * c; // because c = 0.5 * RSS of the (unweighted) residuals
		scale = rss / double(dof);
	} else {
		scale = 1.0; // dof=0
	}
	
	Mat2 Cov = scale * Hinv;

	out.mu_std = std::sqrt(std::max(0.0, Cov(0,0)));
	const double var_s = std::max(0.0, Cov(1,1));
	out.sigma_std = out.sigma * std::sqrt(var_s);

	// a0 uncertainty by propagation: Var(a0) ≈ [da0_dmu da0_ds] * Cov * [da0_dmu da0_ds]ᵀ
	const Vec2 grada0(amp.da0_dmu, amp.da0_ds);
	const double var_a0 = std::max(0.0, (grada0.transpose() * Cov * grada0).value());
	out.a0_std = std::sqrt(var_a0);

	/* Final check, just in case. */
	out.ok = std::isfinite(out.a0) && std::isfinite(out.mu) && std::isfinite(out.sigma) &&
		std::isfinite(out.a0_std) && std::isfinite(out.mu_std) && std::isfinite(out.sigma_std) &&
		(out.sigma > 0.0);

	return out;
}

template<size_t N>
FGResult FastGaussFit (
	const Eigen::Matrix<double, N, 1>& x,
    const Eigen::Matrix<double, N, 1>& y,
    int max_iters   = FG_DEFAULT_MAX_ITERS,
    double y_sigma  = FG_DEFAULT_Y_SIGMA,  // Optionally known noise. for FOOT: around ~2.8 ADC units.
    double lambda0  = FG_DEFAULT_LAMBDA0,  // Lambda parameter in LM algorithm.
    double tol_step = FG_DEFAULT_TOL_STEP,
    double tol_grad = FG_DEFAULT_TOL_GRAD 
) noexcept {
	/* Here we construct initial parameters based on the analytical best guess. */
	Eigen::Matrix<double, 2, 1> init;
	const double ysum = y.sum();
	const double mu = x.dot(y) / ysum; // CoG - centre of gravity
	const auto xc = x.array() - mu;

	const double sigma = sqrt( (y.array() * xc.square()).sum() / ysum ); // ≈ sqrt(Var[CoG])
	init[0] = mu; init[1] = sigma;

#ifdef FGF_LIB_DEBUG_
	std::cout << "Trying with parameters (μ,σ) = " << init.transpose() << std::endl;
#endif

	return FastGaussFit<N>(x, y,
		init,
		max_iters, y_sigma, lambda0, tol_step, tol_grad
	);
}
