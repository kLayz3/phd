#pragma once

#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TF1.h"
#include "Verbosity.hxx"

static constexpr double GAUSS_FIT_SIDE_RATIO_DEFAULT = 1.5;
/**
 * @brief For a projection, TH1D, quickly fit a gauss around it's peak value, and +-2 sigma 
 * @note Thread-unsafe.
 */
inline std::pair <
	std::array<double, 3>,
	std::array<double, 3>
> GaussFitMax (
	TH1D* h, 
	double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	Verbosity v = Verbosity::SILENT
) {
	static uint64_t incrementer_ = 0;

	const double m = h->GetXaxis()->GetBinCenter( h->GetMaximumBin() );
	const double s = h->GetStdDev();
	
	const double fitMin = m - side_ratio*s;
	const double fitMax = m + side_ratio*s;
	if(v > 1)
		printf("[GaussFitMax] (%s) Performing gaus fit around %.2f: [%.2f, %2.f]\n",
			h->GetTitle(), m, fitMin, fitMax);
	TF1 f(Form("f_%lu_%.1f_%.1f", incrementer_++, fitMin, fitMax),"gaus", fitMin, fitMax);
	TFitResultPtr res = h->Fit(&f, "Q0SR");

	if(!res.Get() || (int)res != 0) { // Fit failed.
		if(v > Verbosity::SILENT) {
			fprintf(stderr, "No fit performed...\n");
			fprintf(stderr, "Hist name: '%s', m=%.2f, s=%.2f\n",
				h->GetName(), m, s);
		}
		/* In this case, don't quietly return the NAN's, try to give the best estimate if it
		 * were just a random distribution (uniform) */
		// return {{NAN, NAN, NAN}, {NAN, NAN, NAN}};
		return {{ NAN, m, s }, { NAN, s, s/3 }};
	}
	if(v > 1) std::cout << "[GaussFitMax] (" << h->GetTitle() << ") Result: {" 
		<< f.GetParameter(0) << " +- " << f.GetParError(0) << ", "
		<< f.GetParameter(1) << " +- " << f.GetParError(1) << ", "
		<< f.GetParameter(2) << " +- " << f.GetParError(2) << "}" << std::endl;

	return { 
		std::array<double, 3> { 
			f.GetParameter(0), /* Amplitude */
			f.GetParameter(1), /* Mean */
			f.GetParameter(2)  /* Sigma */
		}, 
		std::array<double, 3> { 
			f.GetParError(0), /* Amplitude */
			f.GetParError(1), /* Mean */
			f.GetParError(2)  /* Sigma */
		}
	};
}

enum class fit_info { PROFILE_MAX, GAUSS_MAX };

[[ maybe_unused ]] static int gCol_ = kRed + 1; /* For Gaussian profile. */
[[ maybe_unused ]] static int pCol_ = kMagenta; /* For standard TProfile. */
