#include "FFT.h"

#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH2D.h"
#include <iomanip>

static const char* module_name_ = "FFTX";

std::ostream& operator<<(std::ostream& os, FFTW::Coeff coeff) {
	os << '[' << coeff.a << ", " << coeff.b << ']';
	return os;
}

void FFTW::ComputeCoeffs(Verbosity v, int Q) {
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
	if(v > 0) this->Dump(Q);

	is_computed = true;
}

void FFTW::Dump(int Q) const {
	std::string s; 
	Q = std::min(Q, (int)coeff.size()); 
	
	auto& os = std::cout;
	os << "[FFTW::Dump(int Q=" << Q << "). N =" << n_full << std::endl << std::setprecision(6);
	os << '[';
	if(Q > 0) os << coeff[0];
	for(int i = 1; i < Q; ++i) {
		os << ", " << coeff[i];
	}
	os << ']' << std::endl;
}

double FFTW::Evaluate(const Foreign& c, const int N, double x, double xmin, double xmax, int Q) {
	const double P = xmax - xmin;
	if(P == 0.0)
		throw std::invalid_argument("static FFTW::Evaluate(): xmax == xmin");

	const int qmax = (int)c.size() - 1; // Max order of the fourier sum
	Q = std::min(Q, qmax);
	if(Q < 0) return NAN;

	const double dx = P / N;
	const double theta = 2.0 * pi * (x - xmin - 0.5*dx ) / P;

	double sum = c[0][0]; // DC component.

	for(int k = 1; k <= Q; ++k) {
		sum += c[k][0] * std::cos(k * theta)
			+  c[k][1] * std::sin(k * theta);
	}
	return sum;
}

double FFTW::Evaluate(double x, double xmin, double xmax, int Q) {
	if (!is_computed)
		ComputeCoeffs(Verbosity::SILENT);

	const double P = xmax - xmin;
	if(P == 0.0)
		throw std::invalid_argument("FFTW::Evaluate(): xmax == xmin");

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

template<fit_info sliceType>
std::tuple <
	FFTW,
	TGraphErrors*,
	TGraph*
> DoFFTW (
	TH2D* h2,
	double x_lo,
	double x_hi,
	double side_ratio,
	int Q,
	Verbosity v
) {
	std::unique_ptr<TH1D> pfx;
	TAxis* xax = h2->GetXaxis();
	x_lo = std::max( xax->GetXmin(), x_lo );
	x_hi = std::min( xax->GetXmax(), x_hi );
	xax->SetRangeUser(x_lo, x_hi);

	if constexpr(sliceType == fit_info::PROFILE_MAX) {
		pfx = std::unique_ptr<TProfile>( h2->ProfileX("pfx__") );
		if(!pfx)
			throw std::runtime_error("DoFFTW(): ProfileX() failed");

		pfx->SetDirectory(nullptr);
	}
	else {
		int ibin_first = xax->GetFirst();
		int ibin_last = xax->GetLast();
		double xfirst = xax->GetBinLowEdge(ibin_first);
		double xlast = xax->GetBinUpEdge(ibin_last);
		int nbins = ibin_last - ibin_first + 1;
		
		if(v > 1) { 
			printf("[DoFFT] %s: from \'%s\'; taking bin indices: [%d,%d], "
				"with lo = %.1f, hi = %.1f, N=%d\n", 
				module_name_, h2->GetTitle(), ibin_first, ibin_last, xfirst, xlast, nbins);
		}

		pfx = std::make_unique<TH1D>("pfx__", "pfx__", 
			nbins, xfirst, xlast);
		pfx->SetDirectory(nullptr);
		for(int ix = ibin_first, i=1; ix <= ibin_last; ++ix, ++i) {
			auto py = std::unique_ptr<TH1D>(h2->ProjectionY(Form("_py__%d", ix), ix, ix));
			py->SetDirectory(nullptr);
			auto [result, err] = GaussFitMax(py.get(), side_ratio, v);
			pfx->SetBinContent(i, result[1]);
			pfx->SetBinError(i, err[1]);
			if(v > 1) {
				printf("[%d]: %.2f +- %.2f\n", i, result[1], err[1]);
			}
		}
	}
	const int N = pfx->GetNbinsX();
	
	TGraphErrors* gerr = new TGraphErrors( pfx.get() ); 
	gerr->SetMarkerStyle(20);
	gerr->SetMarkerSize(1.2);
		
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
		throw std::runtime_error("DoFFTW(): fftw_plan_dft_r2c_1d() failed");
	}

	fftw_execute(plan);
	fftw_destroy_plan(plan);
	fftw_free(in);

	out.ComputeCoeffs(v, Q);
	const int Npts = 70;
	TGraph* g = new TGraph(Npts);
	for(int i=0; i<Npts; ++i) {
		double x0 = x_lo + (i+0.5) * (x_hi- x_lo) / Npts; // centre
		double y0 = out.Evaluate(x0, x_lo, x_hi, Q);
		g->SetPoint(i, x0, y0);
	}
	if constexpr(sliceType == fit_info::PROFILE_MAX) {
		g->SetLineColor(pCol_);
	} else { 
		g->SetLineColor(gCol_);
	}

	g->SetLineWidth(4);

	return { std::move(out), gerr, g };
	// In user code, FFTW instance can directly calculate the FFT'ed value with
	// `.Evaluate(...)` call
}

/* Directly instantiate the templates. 
 * I don't really know if Cling can digest the raw fftw calls.
 * Rather blob everything up here and tuck it in, so that cling doesn't whine. */
template std::tuple <
	FFTW,
	TGraphErrors*,
	TGraph*
>
DoFFTW<fit_info::PROFILE_MAX>(TH2D*, double, double, double, int, Verbosity);

template std::tuple <
	FFTW,
	TGraphErrors*,
	TGraph*
>
DoFFTW<fit_info::GAUSS_MAX>(TH2D*, double, double, double, int, Verbosity);
