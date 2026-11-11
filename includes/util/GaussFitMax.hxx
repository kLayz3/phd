#pragma once

/* Explicitly avoiding all monad dependency... */
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TF1.h"
#include "Verbosity.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <utility>

constexpr double GAUSS_FIT_SIDE_RATIO_DEFAULT = 1.5;

using GaussFitMaxResults = std::pair<std::array<double, 3>, std::array<double, 3>>;

/* Fit a Gaussian around the dominant peak in the histogram's selected axis range.
 * The initial width uses the selected peak's local FWHM, with histogram stddev as fallback.
 * niter is a maximum: each iteration uses the previous fitted mean and sigma for its range.
 * Returns { {amplitude, mean, sigma}, {their standard errors} }, or all NaNs on failure.
 * With likelihood=true, use Poisson likelihood for unweighted, nonnegative count histograms.
 * Thread-unsafe ROOT fitting API. Kinda slow function, but very hardened... */
inline GaussFitMaxResults GaussFitMax (
    TH1D* h,
    double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
    uint32_t niter = 2,
    Verbosity v = Verbosity::SILENT,
    bool likelihood = false
) {
	constexpr const char* color_red = "\e[1;91m";
	constexpr const char* color_grn = "\e[1;92m";
	constexpr const char* color_cyn = "\e[1;96m";
	constexpr const char* color_off = "\e[0m";

    constexpr double invalid = std::numeric_limits<double>::quiet_NaN();

    const auto failure = [&](const char* reason) {
        if(v > Verbosity::SILENT)
            std::fprintf(stderr, "%s%s%s:%s%d%s failure in histogram %s%s%s : %s\n",
				color_grn,  __FILE_NAME__, color_off, color_cyn, __LINE__, color_off,
				color_red, h? h->GetName(): "<null>", color_off, reason);
        return std::make_pair(std::array<double, 3>{invalid, invalid, invalid},
                              std::array<double, 3>{invalid, invalid, invalid});
    };
    /* [1] Sanity check */
    if(!h || !std::isfinite(side_ratio) || side_ratio <= 0 || niter == 0)
        return failure("invalid histogram, side ratio, or iteration count");

    const auto* axis = h->GetXaxis();
    const int first = std::max(1, axis->GetFirst());
    const int last = std::min(h->GetNbinsX(), axis->GetLast());
    if(last - first + 1 < 4)
        return failure("fewer than four selected bins");

    for(int bin = first; bin <= last; ++bin) {
        const double content = h->GetBinContent(bin), error = h->GetBinError(bin);
        if(!std::isfinite(content) || !std::isfinite(error) || error < 0 || (likelihood && content < 0))
            return failure("invalid bin contents or uncertainties");
    }

    // [2] Find the largest three-bin weighted average: weights {1,2,1} / 4.
    // Neighboring bins support a real peak, while an isolated spike gets washed out.
    // At the selected range edges, renormalize using only available neighbors.
    // This smoothing is only for initialization; ROOT fits the original histogram.
    const auto smooth = [&](int bin) {
        double sum = 2 * h->GetBinContent(bin), weight = 2;
        if(bin > first) { sum += h->GetBinContent(bin - 1); ++weight; }
        if(bin < last) { sum += h->GetBinContent(bin + 1); ++weight; }
        return sum / weight;
    };
    int peak = first;
    double height = smooth(first), baseline = height;
    for(int bin = first + 1; bin <= last; ++bin) {
        const double value = smooth(bin);
        baseline = std::min(baseline, value);
        if(value > height) { peak = bin; height = value; }
    }
    if(height <= 0 || height <= baseline)
        return failure("no resolved peak");

    // [3] Seed the mean at that bin and the amplitude from its original content.
    // Estimate the baseline as the minimum smoothed content in the selected range.
    // Walk outward to the first crossings of half the peak height above this baseline.
    // Interpolate each crossing between bin centers, then use sigma = FWHM / 2.35482.
    // If either crossing is missing (e.g. a clipped peak), retain GetStdDev() instead.
    double m = axis->GetBinCenter(peak);
    double a = h->GetBinContent(peak);
    double s = h->GetStdDev();
    const double half = baseline + 0.5 * (height - baseline);
    int left = peak, right = peak;
    while(left > first && smooth(left) > half) --left;
    while(right < last && smooth(right) > half) ++right;
    if(left < peak && right > peak && smooth(left) <= half && smooth(right) <= half) {
        const auto crossing = [&](int outside, int inside) {
            const double y0 = smooth(outside), y1 = smooth(inside);
            return axis->GetBinCenter(outside) + (half - y0) / (y1 - y0)
                * (axis->GetBinCenter(inside) - axis->GetBinCenter(outside));
        };
        s = (crossing(right, right - 1) - crossing(left, left + 1)) / 2.3548200450309493;
    }
    if(!std::isfinite(s) || s <= 0 || !std::isfinite(a) || a <= 0)
        return failure("cannot estimate a positive peak width and amplitude");

    const double bin_width = axis->GetBinWidth(peak);
    const double lower = axis->GetBinLowEdge(first), upper = axis->GetBinUpEdge(last);
    std::pair<std::array<double, 3>, std::array<double, 3>> result;

	// [4] Refit the core inside mean +/- side_ratio*sigma, clipped to the selected range.
    // Each successful fit supplies the parameters and window for the next iteration.
    for(uint32_t iter = 0; iter < niter; ++iter) {
        const double fit_min = std::max(lower, m - side_ratio * s);
        const double fit_max = std::min(upper, m + side_ratio * s);
        // Three free parameters need more than three usable bins. Also require four
        // positive-content bins so empty bins alone cannot make a noise spike look resolved.
        // Chi-square excludes zero-error bins; likelihood includes empty count bins.
        int usable = 0, populated = 0;
        for(int bin = first; bin <= last; ++bin) {
            const double x = axis->GetBinCenter(bin);
            if(x < fit_min || x > fit_max) continue;
            if(likelihood || h->GetBinError(bin) > 0) ++usable;
            if(h->GetBinContent(bin) > 0) ++populated;
        }
        if(usable < 4 || populated < 4) {
            if(iter > 0) break; // Keep the last valid fit rather than fit an unresolved window.
            return failure("peak window has fewer than four usable, populated bins");
        }
        if(v > Verbosity::INFO)
            std::fprintf(stderr, "[GaussFitMax: %u/%u] %s: range [%.6g, %.6g], seed {%.6g, %.6g, %.6g}\n",
                         iter + 1, niter, h->GetName(), fit_min, fit_max, a, m, s);

        // [5] Explicitly reuse the seed and constrain the fit to a positive-width peak
        // centered inside this window. The amplitude cap is a kinda runaway guard.
        // Q=quiet, N=no storage/drawing, S=return fit result, R=use function range,
        // B=keep our initialization/limits for predefined gaus; L=Poisson likelihood.
        TF1 f("gauss_fit_max", "gaus", fit_min, fit_max, TF1::EAddToList::kNo);
        f.SetParameters(a, m, s);
        f.SetParLimits(0, 0, 100 * height);
        f.SetParLimits(1, fit_min, fit_max);
        f.SetParLimits(2, 0.1 * bin_width, upper - lower);
        const TFitResultPtr fit = h->Fit(&f, likelihood ? "QNSRBL" : "QNSRB");

		// [6] A zero minimizer status alone is insufficient: require a valid result
        if(!fit.Get() || int(fit) != 0 || !fit->IsValid() || fit->CovMatrixStatus() < 2 || fit->Ndf() <= 0)
            return failure("fit failed or covariance is not usable?");

        for(int i = 0; i < 3; ++i) {
            result.first[i] = f.GetParameter(i);
            result.second[i] = f.GetParError(i);
            if(!std::isfinite(result.first[i]) || !std::isfinite(result.second[i]) || result.second[i] < 0)
                return failure("non-finite fit parameters or uncertainties");
        }
        // Reject nonpositive amplitudes, means pinned to a range edge, and widths
        // below half the seed bin width or as broad as the entire selected range.
        const auto& p = result.first;
        if(p[0] <= 0 || p[1] <= fit_min || p[1] >= fit_max || p[2] < 0.5 * bin_width || p[2] >= upper - lower)
            return failure("nonphysical, unresolved, or boundary-limited fit");
        if(v > Verbosity::INFO)
            std::fprintf(stderr, "[GaussFitMax: %u/%u] %s: {%.6g, %.6g, %.6g}, errors {%.6g, %.6g, %.6g}\n",
                         iter + 1, niter, h->GetName(), p[0], p[1], p[2],
                         result.second[0], result.second[1], result.second[2]);

        // [7] Stop when the mean changes by <= 1% of a seed bin and sigma changes by
        // <= 0.1% of max(previous sigma, seed bin width). Otherwise repeat, up to niter.
        const bool stable = std::abs(p[1] - m) <= 0.01 * bin_width
                         && std::abs(p[2] - s) <= 0.001 * std::max(s, bin_width);
        a = p[0]; m = p[1]; s = p[2];
        if(stable) break;
    }
    // Full covariant matrix I don't really care about. Trust the gauss-chan 🥺 👉👈
    return result;
}

inline bool VerifyGaussFitMax(const GaussFitMaxResults& fit) noexcept {
	return std::isfinite(fit.first[0]) &&
		std::isfinite(fit.first[1]) &&
		std::isfinite(fit.first[2]) &&
		std::isfinite(fit.second[0]) &&
		std::isfinite(fit.second[1]) &&
		std::isfinite(fit.second[2]) &&
		fit.first[0] > 0.0 &&
		fit.first[2] > 0.0;
}

enum class fit_info { PROFILE_MAX, GAUSS_MAX };

[[ maybe_unused ]] static int gCol_ = kRed + 1; /* For Gaussian profile. */
[[ maybe_unused ]] static int pCol_ = kMagenta; /* For standard TProfile. */
