#include "FFT.h"

#include "TProfile.h"
#include "TH2D.h"

void FFTW::ComputeCoeffs() {
	// k = 0
	coeff[0].a = data[0][0] / n_full;
	coeff[0].b = 0.0;
	// 1 <= k < n_full/2
	for(int k = 1; k < n_freq; ++k) {
		const double Re = data[k][0];
		const double Im = data[k][1];

		// Nyquist theorem :-)
		if ((n_full % 2 == 0) && (k == n_full / 2)) {
			coeff[k].a = Re / n_full; 
			coeff[k].b = 0.0;
		} else {
			coeff[k].a =  2.0 * Re / n_full;
			coeff[k].b = -2.0 * Im / n_full;
		}
	}

	is_computed = true;
}
double FFTW::Evaluate(double x, double xmin, double xmax, int Q) {
	if (!is_computed)
		ComputeCoeffs();

	const double P = xmax - xmin;
	if(P == 0.0)
		throw std::invalid_argument("Evaluate(): xmax == xmin");

	const int qmax = n_full / 2;
	Q = std::min(Q, qmax);

	const double dx = P / n_full;
	const double theta = 2.0 * pi * (x - xmin - 0.5*dx ) / P;

	double sum = coeff[0].a;

	for(int k = 1; k <= Q; ++k) {
		sum += coeff[k].a * std::cos(k * theta)
			+ coeff[k].b * std::sin(k * theta);
	}
	return sum;
}

FFTW DoFFTX(TH2D* h2) {
	auto pfx = std::unique_ptr<TProfile>( h2->ProfileX("pfx") );
	if(!pfx) {
		throw std::runtime_error("DoFFTX(): ProfileX() failed");
	}
	pfx->SetDirectory(nullptr);

	const int N = pfx->GetNbinsX();

	double* in = (double*) fftw_malloc(sizeof(double) * N);
	if(!in) {
		throw std::bad_alloc{};
	}
	
	for(int i=1; i <= N; ++i)
		in[i-1] = pfx->GetBinContent(i);

	FFTW out{N}; 

	fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out.data, FFTW_ESTIMATE);
    if(!plan) {
		fftw_free(in);
		throw std::runtime_error("DoFFTX(): fftw_plan_dft_r2c_1d() failed");
	}

	fftw_execute(plan);
	fftw_destroy_plan(plan);
	fftw_free(in);

	out.ComputeCoeffs();
	return out;
	// In user code, FFTW instance can directly calculate the FFT'ed value with
	// `.Evaluate(...)` call
}

/* To check if it works fine:

const double xmin = pfx->GetXaxis()->GetXmin();
const double xmax = pfx->GetXaxis()->GetXmax();
const int Q = 4;

TF1* fFourier = new TF1(
    "fFourier",
    [&fft, xmin, xmax, Q](double* xx, double*) {
        return fft.Evaluate(xx[0], xmin, xmax, Q);
    },
    xmin, xmax, 0
);

fFourier->SetLineColor(kRed);
fFourier->SetLineWidth(3);

h2->Draw("COLZ");
pfx->SetMarkerStyle(20);
pfx->SetMarkerSize(0.7);
pfx->SetMarkerColor(kBlack);
pfx->Draw("SAME");

fFourier->Draw("SAME");

*/
