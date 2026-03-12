#pragma once

#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TF1.h"

/**
 * @brief For a projection, TH1D, quickly fit a gauss around it's peak value, and +-2 sigma 
 * @note Thread-unsafe.
 */
namespace mnd {
inline
std::array<double, 3> GaussFitMax(TH1D* h) {
	double m = h->GetXaxis()->GetBinCenter( h->GetMaximumBin() );
	double s = h->GetStdDev();
	
	double fitMin = m - 2*s;
	double fitMax = m + 2*s;
	TF1 f("f","gaus", fitMin, fitMax);
	TFitResultPtr res = h->Fit(&f, "Q0S");

	if(!res.Get() || (int)res != 0) { // Fit failed.
		fprintf(stderr, "No fit performed...\n");
		fprintf(stderr, "Hist name: '%s', m=%.2f, s=%.2f\n",
			h->GetName(), m, s);
		return {NAN, NAN, NAN};
	}
	
	return { 
		f.GetParameter(0),
		f.GetParameter(1),
		f.GetParameter(2)
	};
}
}
