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
	static constexpr double pi = 3.141592653589793238462643383279502884;
	struct Coeff { double a, b; };

	int n_full{0};   // number of real input samples
	int n_freq{0};   // number of r2c output bins = n_full/2 + 1
	
	fftw_complex *data {}; 
	std::vector<Coeff> coeff {};
	bool is_computed {false};

	FFTW(int n_samples) : 
		n_full(n_samples),
		n_freq(n_samples / 2 + 1), // Nyquist theorem :-)
		data( (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* n_freq) ),
		coeff(n_freq)
	{
		if(!data) 
			throw std::bad_alloc{};
	}

	~FFTW() {
		if(data) fftw_free(data);
	}
	
	FFTW(const FFTW&) = delete;
	FFTW& operator=(const FFTW&) = delete;

	FFTW(FFTW&& rhs) noexcept : 
		n_full(rhs.n_full),
		n_freq(rhs.n_freq),
		data(rhs.data),
		coeff(std::move(rhs.coeff)),
		is_computed(rhs.is_computed)
	{
		rhs.n_full = 0;
		rhs.n_freq = 0;
		rhs.data = nullptr;
		rhs.is_computed = false;
	}

	FFTW& operator=(FFTW&& rhs) noexcept {
		if(this != &rhs) {
			if(data) fftw_free(data);

			n_full = rhs.n_full;
			n_freq = rhs.n_freq;
			data = rhs.data;
			coeff = std::move(rhs.coeff);
			is_computed = rhs.is_computed;

			rhs.n_full = 0;
			rhs.n_freq = 0;
			rhs.data = nullptr;
			rhs.is_computed = false;
		}
		return *this;
	}
	
	void ComputeCoeffs(Verbosity = Verbosity::SILENT, int = Q_DEFAULT);
	void Dump(int = Q_DEFAULT) const;

	// Evaluate in the sample-aligned basis:
	// x_j = xmin + (j + 1/2) * (P/n_full)
	double Evaluate(double x, double xmin, double xmax, int Q = INT_MAX);
	
	/* A static function to evaluate based on `std::vector<std::array<double, 2>>` packed- 
	 * foreign coefficients, and `N` supplied. */
	using Foreign = std::vector<std::array<double, 2>>;
	static double Evaluate(const Foreign& c, const int N, double x, double xmin, double xmax, int Q = INT_MAX);

	inline double DC() const noexcept {
		return coeff[0].a;	
	}
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
