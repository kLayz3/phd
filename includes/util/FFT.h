#pragma once

/* Make sure to link against libfftw3.so,
 * it's not part of libc! */

#include <climits>
#include <new>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>

#include "GaussFitMax.hxx"
#include "Verbosity.hxx"

extern "C" {
	#include <fftw3.h>
}

/* Make it RAII friendly.
 * @note Thread un-safe. */
struct FFTW {
	static constexpr int Q_DEFAULT = 10;
	static constexpr double pi = M_PI;
	struct Coeff { double a, b; };

	int n_full{0};   // number of real input samples
	int n_freq{0};   // number of r2c output bins = n_full/2 + 1
	
	fftw_complex *data = nullptr;
	std::vector<Coeff> coeff {};
	bool is_computed {false};

	FFTW(int );
	
	FFTW(const FFTW&) = delete;
	FFTW& operator=(const FFTW&) = delete;

	FFTW(FFTW&& rhs) noexcept;
	FFTW& operator=(FFTW&& rhs) noexcept;
	~FFTW();
	
	void ComputeCoeffs(Verbosity = Verbosity::SILENT, int = Q_DEFAULT);
	void Dump(int = Q_DEFAULT) const;

	// Evaluate in the sample-aligned basis:
	// x_j = xmin + (j + 1/2) * (P/n_full)
	double Evaluate(double x, double xmin, double xmax, int Q = INT_MAX);
	
	/* A static function to evaluate based on `std::vector<std::array<double, 2>>` packed
	 * foreign coefficients, and `N` supplied. */
	using Foreign = std::vector<std::array<double, 2>>;
	static double Evaluate(const Foreign& c, const int N, double x, double xmin, double xmax, int Q = INT_MAX);

	double DC() const noexcept;
};

class TH2D;
class TGraph;
class TGraphErrors;

template<fit_info>
std::tuple <
	FFTW,
	TGraphErrors*,
	TGraph*
> DoFFTW (
	TH2D*,
	double = -DBL_MAX,
	double = DBL_MAX,
	double = GAUSS_FIT_SIDE_RATIO_DEFAULT,
	int = 6,
	Verbosity = Verbosity::SILENT
);

std::ostream& operator<<(std::ostream& , const FFTW::Coeff& );
