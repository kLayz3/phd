#include "TFOOTMapProc.h"
#include "TFRSGo4Cont.hxx"
#include "TF1.h"
#include "TGraph.h"
#include "eigen/Eigen/Core"
#include "libs.hh"
#include <cassert>
#include <cfloat>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <cmath>
#include "AuxFunctions.hh"

#include "eigen/Eigen/Dense"

#if defined(__GNUC__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Winvalid-offsetof"
#endif

static_assert(TFOOTMapProc::N_STRIPS == TFOOTMapProc::N_ASIC * TFOOTMapProc::N_STRIPS_PER_ASIC, "Failed build: nstrip != nasic*nstrip_per_asic!\n");

using nlohmann::json;
json TFOOTMapProc::_bad_strips{};

	//struct FOOTView {
	//	bool* TSBad; u32* TLO; u32* THI; u32* SY;
	//	u32* Ndata; u32* E;
const TFOOTMapProc::FOOTView TFOOTMapProc::GetPtrs(const TFRSSortEvent* e, int N) {
#define GET_FOOT_GO4(N) \
	case(N): { \
		return { (bool)e->FOOT##N##TSBAD, \
			     (u32 )e->FOOT##N##TLO  , \
			     (u32 )e->FOOT##N##THI  , \
			     (u32 )e->FOOT##N##SY   , \
			     (u32 )e->FOOT##N       , \
			     (u32*)e->FOOT##N##E   }; \
		}
	switch(N) {
		GET_FOOT_GO4( 1)
		GET_FOOT_GO4( 2)
		GET_FOOT_GO4( 3)
		GET_FOOT_GO4( 4)
		GET_FOOT_GO4( 5)
		GET_FOOT_GO4( 6)
		GET_FOOT_GO4( 7)
		GET_FOOT_GO4( 8)
		GET_FOOT_GO4( 9)
		GET_FOOT_GO4(10)
		GET_FOOT_GO4(11)
		GET_FOOT_GO4(12)
		GET_FOOT_GO4(13)
		GET_FOOT_GO4(14)
		GET_FOOT_GO4(15)
		GET_FOOT_GO4(16)
		GET_FOOT_GO4(17)
		GET_FOOT_GO4(18)
		GET_FOOT_GO4(19)
		GET_FOOT_GO4(20)
		GET_FOOT_GO4(21)
		GET_FOOT_GO4(22)
		GET_FOOT_GO4(23)
		GET_FOOT_GO4(24)
		GET_FOOT_GO4(25)
		GET_FOOT_GO4(26)
		GET_FOOT_GO4(27)
		GET_FOOT_GO4(28)
		GET_FOOT_GO4(29)
		
		default:
			ERROR("Request for member %d, not in [1,29]\n", N);
	}
}

/* Assumption is that each bin is exactly 1 units wide. Can work also generally,
 * but has to be tweaked a bit. */

TFOOTMapProc::GaussFitParams TFOOTMapProc::FitGauss(const TH1D* h) {
	constexpr double RATIO_THR = 0.2;

	/* First find the histogram max bin. */
	const int maxBin = h->GetMaximumBin();
	const double binWidth2 = h->GetBinWidth(maxBin) / 2;

	const double N0 = h->GetBinContent(maxBin);

	/* Idea is to fit a parabola to log(count), around 1-2 sigma. */
	std::vector<std::pair<double, double>> points;
	points.emplace_back(maxBin + binWidth2, std::log(N0));

	int i;
	
	i= maxBin + 1;
	while( h->GetBinContent(i) / N0 > RATIO_THR ) {
		points.emplace_back(i + binWidth2, std::log(h->GetBinContent(i)) );
		++i;
	}
	i = maxBin - 1;
	while( h->GetBinContent(i) / N0 > RATIO_THR ) {
		points.emplace_back(i + binWidth2, std::log(h->GetBinContent(i)) );
		--i;
	}

	auto [a,b,c] = FitParabolaLeastSquares(std::move(points));
	
	return { - b / (2*c) , 1 / sqrt(-2*c) }; 
}

// To fit `n` points (xi,yi) with parabola, solve the following 3x3:
/* ( n      Σxi    Σxi^2 ) ( a )   ( Σyi     )
 * ( Σxi    Σxi^2  Σxi^3 ) ( b ) = ( Σxiyi   )
 * ( Σxi^2  Σxi^3  Σxi^4 ) ( c )   ( Σxi^2yi )
 */
TFOOTMapProc::ParabolaFitParams
TFOOTMapProc::FitParabolaLeastSquares(TFOOTMapProc::Points&& points) {
	if(points.size() < 3) return { std::nan(""), std::nan(""), std::nan("") };

	double Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0;
	double Sy = 0, Sxy = 0, Sx2y = 0;
	const int n = (int)points.size();
	
	for(int i=0; i<n; ++i) {
		auto [x,y] = points[i];
		double x2 = x * x;
		Sx += x;
		Sx2 += x2;
		Sx3 += x2 * x;
		Sx4 += x2 * x2;
		Sy += y;
		Sxy += x * y;
		Sx2y += x2 * y;
	}
	Eigen::Matrix3d M;
	M << n,   Sx,  Sx2,
	     Sx,  Sx2, Sx3,
	     Sx2, Sx3, Sx4;
	
	Eigen::Vector3d v(Sy, Sxy, Sx2y);

	M.colPivHouseholderQr().solve(v);
	return { v(0,0), v(1,0), v(2,0) };
};

void TFOOTMapProc::ProcessEntry() noexcept {
	switch(process_type) {
		case kINITIAL_BATCH: 
			ProcessInitialPedestal();
			break;
		case kMAIN:
			ProcessEventPedestal();
			break;
	}
}

void TFOOTMapProc::CalcGlobalPedestal() {
	TH2I* h = out.h2_raw_tmp;
	if(h->GetEntries() == 0) {
		ERROR("Ran over the TTree initial batch, but found 0 events with data? " EMPH(FOOT: %d) ", Setting all pedestals to 0.", N);
		return;
	}
	FOR(i, N_STRIPS) {
		auto slice = std::unique_ptr<TH1D>( h->ProjectionY(_MSG("_py%d-%d", N, i), i+1, i+1) );
		slice->SetDirectory(nullptr);
		
		auto [pedestal, sigma0] = FitGauss(slice.get());
		out.h2_ped0->Fill(i, pedestal);
		out.h2_sigma0->Fill(i, sigma0);
		
		current_gped[i] = pedestal;
	}
	h->Reset("ICESM");
}

/* Ok weirdest hack ever, how to get pointer access to .FOOT##N##E if `N`
 * isn't in a preproc or hardcoded? Behold, breaking the universe, and my own sanity.
 * (FYI I wrote the unpacker, too ¯\_(ツ)_/¯ ) */
void TFOOTMapProc::ProcessInitialPedestal() noexcept {
	if(N <= 0) ERROR("Starting to process, but local foot label <= 0 ?");

	const TFRSSortEvent* sortev = std::get<0>( this->in ).raw();
	
	auto [_0, _1, _2, _3, footN, data] = GetPtrs(sortev, N);
	if(footN == 0) return; 
	
	int iraw;	
	FOR(i, N_STRIPS) {
		iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;
		out.h2_raw_tmp->Fill(i, data[iraw]);
	}
}

//#define CALC_OFFSET_FROM_MEDIAN

void TFOOTMapProc::ProcessEventPedestal() noexcept {
	TFRSGo4Cont & input  = std::get<0>( this->in );
	out.Clean();

	const TFRSSortEvent* sortev = &input.inner();

	auto [_0, _1, _2, _3, footN, data] = GetPtrs(sortev, N);
	if(footN == 0) return;
		

	/* Subtract the global pedestal. */
	/* Algorithm (per groups of 64-strips) is the following: j=0,1,2, ... 63
	 * (1) Calculate the offsets of each: _FOOTE[j] - gped[j] and store in a 64-length array.
	 * (2) Sort the array.
	 * (3) Trim the lowest and highest N_TRIM number of entries
	 * (4) From what is remaining, (re)calculate the mean, yielding the offset `o[asic]` of each asic
	 * (5) Offset the calculation yielding the final calibrated ADC value:
	 * --> adc[i] = _FOOTE[i] - gped[i] - o[i / N_STRIPS_PER_ASIC] */

	/* Should be between 5..10 */
#define N_TRIM_FINE_PED_LO 8
#define N_TRIM_FINE_PED_HI 8

	double ped_off_med, ped_off_avg, adc_final; 
	
	FOR(asic, N_ASIC) /* 0..=9 */ {
		const int i0 = asic * N_STRIPS_PER_ASIC;
		int i, iraw; 
		
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip; // <--- 'true strip'
			iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;
			
			out.h2_raw_tmp->Fill(i, data[iraw]);
			
			/* Calculate the systematic shift per asic (group of 64 strips). */
			ped_offset[strip] = data[iraw] - current_gped[i];
		}

		/* Initial strip #0 is uncoupled. It contributes nothing. 
		 * Just remove it from the analysis. */
		ped_offset[0] = -DBL_MAX;

		std::sort(ped_offset.begin(), ped_offset.end());

		ped_off_med = util::median( ped_offset ); 

		ped_off_avg = std::accumulate (
			ped_offset.cbegin() + N_TRIM_FINE_PED_LO, 
			ped_offset.cend()   - N_TRIM_FINE_PED_HI, (double)0.0
		);
		ped_off_avg /= (N_STRIPS_PER_ASIC - N_TRIM_FINE_PED_LO - N_TRIM_FINE_PED_HI);

		out.h2_ped_off_med->Fill(asic, ped_off_med);
		out.h2_ped_off_avg->Fill(asic, ped_off_avg);
		out.h2_ped_off_diff->Fill(asic, ped_off_avg - ped_off_med);
	
		/* Repeat the loop, subtract the collective small offset (fine correction). */
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip;
			int iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;

			adc_final = data[iraw] - current_gped[iraw];
			
			if(strip > 0) {
				adc_final -=
#ifdef CALC_OFFSET_FROM_MEDIAN
				ped_off_med
#else
				ped_off_avg
#endif
				; /* The initial first strip in ASIC is uncoupled, just let it be. */
			} 
			else { /* Just for the initial strip that's uncoupled, wash away binning all the values into a single bin. */
				adc_final += rand() / (double)RAND_MAX ;
			}
			out.h2_corr->Fill(i, adc_final);
			out.inner().FOOTE[i] = adc_final;
		}
	}

	/* If we sample enough events, recalculate the global pedestal quickly from the batch. */
	if((++nsampled) == N_BATCH_PEDESTAL) {
		CalcGlobalPedestal();	
		nsampled = 0;
	}
}

int TFOOTMapProc::ParseStaticBadStrips() {
	TFOOTMapCont& output = this->out;
	if(! output.bad_strips)
		ERROR("Output bad strips container uninitialized (FOOT_N=%d). Did you call TFOOTMapCont::Init(..) beforehand?",
			N);

	const char* key = Form("FOOT%d", N);
	std::vector<int> parsed{};
	if(auto it = _bad_strips.find(key); it != _bad_strips.end())
		util::parse_json_as_int_vec(parsed, *it);
	
	std::vector<int>& v = *output.bad_strips; 
	v.insert(v.end(),
		parsed.begin(),
		parsed.end()
	);

	/* Remove duplicates and sort it.
	 * Faster than making a hashmap via ctor, 
	 * check: https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
	 * */
	std::unordered_set<int> s;
	for(int i: v) s.insert(i);
	v.assign(s.begin(), s.end());
	std::sort(v.begin(), v.end());

	if(v.size() > 0) {
		WARN("[%s] Found " EMPH(%zu) " bad strips (%zu from param files): ", key, v.size(), parsed.size());
		for(auto x : v) printf("%s%d%s ", BOLD, x, KNRM);
		printf("\n");
	}
	
	return (int)v.size();
}

void TFOOTMapProc::CalcFinalPedestal() {
	if(out.h2_corr->GetEntries() == 0) {
		WARN("Ran over the TTree, but found 0 events with calibrated data?" EMPH(FOOT: %d\n), N);
		return;
	}
	FOR(i, N_STRIPS) {
		auto slice = std::unique_ptr<TH1D>( out.h2_corr->ProjectionY(_MSG("_py%d-%d", N, i), i+1, i+1) );
		slice->SetDirectory(nullptr);
		
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 10);
		double fitMax = slice->GetBinCenter(maxBin + 10);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		
		out.gped_sf->at(i) = fitF->GetParameter(2);
		out.gr_s1->SetPoint(i, i, out.gped_sf->at(i));
	
		if(fitF->GetParameter(2) > BAD_STRIP_CUTOFF_HI or 
			fitF->GetParameter(2) < BAD_STRIP_CUTOFF_LO) {
			out.bad_strips->push_back(i);
		}
	}
	ParseStaticBadStrips();
}


