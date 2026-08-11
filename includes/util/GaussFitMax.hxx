#pragma once

#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TF1.h"
#include "Verbosity.hxx"

static constexpr double GAUSS_FIT_SIDE_RATIO_DEFAULT = 1.5;
/**
 * @brief For a projection, TH1D, quickly fit a gauss around it's peak value, and +-2 sigma.
 * Returns a pair, where first element is the calculated {A,μ,σ} and second element its corresponding variances
 * (diagonal of the full covariant matrix).
 * @note Thread-unsafe. The predicate that fit succeeded is: `std::isfinite( result.first[0] )`
 */
inline std::pair <
	std::array<double, 3>,
	std::array<double, 3>
> GaussFitMax (
	TH1D* h,
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	uint32_t niter = 2,
	Verbosity v = Verbosity::SILENT
) {
	assert(niter > 0 && "Must be at least 1 iteration passed here.");
	static uint64_t incrementer_ = 0;
	double s = NAN, m = NAN, a = NAN;
	double ss, ms, as; 

	for(uint32_t iter=0; iter < niter; ++iter) {
		m = std::isnan(m) ? h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ) : m;
		s = std::isnan(s) ? h->GetStdDev() : s;
		const double fitMin = m - side_ratio*s;
		const double fitMax = m + side_ratio*s;
		if(v > 1)
			fprintf(stderr, "[GaussFitMax: %u/%u] (%s) Performing gaus fit around %.2f: [%.2f, %2.f]\n",
				iter+1, niter, h->GetTitle(), m, fitMin, fitMax);
		TF1 f(Form("f_%lu_%.1f_%.1f", incrementer_++, fitMin, fitMax),"gaus", fitMin, fitMax);
		TFitResultPtr res = h->Fit(&f, "Q0SR");

		if(!res.Get() || (int)res != 0) { // Fit failed.
			if(v > Verbosity::SILENT) {
				fprintf(stderr, "[GaussFitMax: %u/%u] (%s) no fit performed. Log: m=%.2f, s=%.2f\n ",
					iter+1, niter, h->GetTitle(), m, s);
			}
			/* In case the fit diverges, don't quietly return the NAN's,
			 * try to give the best estimate if it were just a random distribution (uniform). */
			return {{ NAN, m, s }, { NAN, s, s/3 }};
		}
		if(v > 1) { 
			fprintf(stderr, "[GaussFitMax: %u/%u] (%s) result: {%.2f ± %.2f, %.2f ± %.2f, %.2f ± %.2f}\n",
				iter+1, niter, h->GetTitle(), 
				f.GetParameter(0), f.GetParError(0),
				f.GetParameter(1), f.GetParError(1),
				f.GetParameter(2), f.GetParError(2));
		}
		a = f.GetParameter(0);
		m = f.GetParameter(1);
		s = f.GetParameter(2);
		as = f.GetParError(0);
		ms = f.GetParError(1);
		ss = f.GetParError(2);
	}

	return { 
		std::array<double, 3> {
			a, /* Amplitude */
			m, /* Mean */
			s  /* Sigma */
		}, 
		std::array<double, 3> {
			as, /* Amplitude */
			ms, /* Mean */
			ss  /* Sigma */
		}
		// Full covariant matrix I don't really care about. Trust the gauss-chan 🥺 👉👈
	};
}

enum class fit_info { PROFILE_MAX, GAUSS_MAX };

[[ maybe_unused ]] static int gCol_ = kRed + 1; /* For Gaussian profile. */
[[ maybe_unused ]] static int pCol_ = kMagenta; /* For standard TProfile. */
