#include "TFOOTMapProc.h"
#include "TFOOTMapCont.h"
#include "TFRSGo4Cont.hxx"

#include "TF1.h"
#include "TGraph.h"
#include "TParameter.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <cmath>
#include <regex>
#include <mutex>

static_assert(TFOOTMapProc::N_STRIPS == TFOOTMapProc::N_ASIC * TFOOTMapProc::N_STRIPS_PER_ASIC, "Failed build: nstrip != nasic*nstrip_per_asic!\n");

static inline void Trim(std::string& s) noexcept {
	s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c){return std::isspace(c);}), s.end());	
}

using nlohmann::json;
json TFOOTMapProc::_bad_strips{};

static std::mutex g_foot_fit_mtx{};

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

TFOOTMapProc::TFOOTMapProc(TFOOTMapCont& out, const TFRSGo4Cont& in, 
	NBatchPedestal n_batch, CableSwapped is_cabling_swapped, double dt_veto_) : Base(out, in), 
	is_swapped(is_cabling_swapped),
	n_batch_pedestal(n_batch),
	dt_veto(1000 * dt_veto_),
	N(out.FOOT_N) 
	{
		/* Write out the `dt_veto` into output TOnce object. */
		out.dt_veto->SetVal(dt_veto_);
	}

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

/* By default, this isn't thread safe, so we guard it. */
void TFOOTMapProc::CalcGlobalPedestal() {
	TH2I* h = out.h2_raw_tmp;
	if(h->GetEntries() == 0) {
		WARN("Ran over the raw Go4 FOOT data batch, but found 0 events with data? " EMPH(FOOT: %d) ", Setting keeping old pedestals / setting to 0.", N);
		return;
	}

	for(int i=0; i < N_STRIPS; ++i) {
		auto slice = std::unique_ptr<TH1D>( h->ProjectionY(mnd::msg("_py%d-%d", N, i), i+1, i+1) );
		slice->SetDirectory(nullptr);
		
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 10);
		double fitMax = slice->GetBinCenter(maxBin + 10);

		if(slice->GetEffectiveEntries() > 3)
			slice->Fit("gaus", "Q", "", fitMin, fitMax);

		TF1* fitF = slice->GetFunction("gaus"); // Quiet null if no fit performed for empty data.
		
		double pedestal, sigma;
		if(fitF == nullptr) {
			//WARN("FOOT%d, strip %d, couldn't fit Gauss for y-projection of raw data? Skipping. "
			//	"Debug: maxBin = %d, integral = %.1f\n", N, i, maxBin, slice->GetEntries());
			pedestal = slice->GetBinContent(maxBin);
			sigma    = BAD_STRIP_CUTOFF_HI + 10;
		} else {
			pedestal = fitF->GetParameter(1);
			sigma    = fitF->GetParameter(2);
		}
	
		out.h2_gped->Fill(i, pedestal);
		out.h2_gped_sigma->Fill(i, sigma);
	
		current_gped[i] = pedestal;
	}

	h->Reset("ICESM");

	if(!initial_calculated) {
		memcpy(initial_gped.data(), current_gped.data(), N_STRIPS * sizeof(double));
		initial_calculated = true;
	}

	for(int i=0; i < N_STRIPS; ++i) {
		out.h2_gped_per_batch->Fill (
			batch_index,
			i * 2 * TFOOTMapCont::N_GPED_CHANGE_TOLERANCE 
			+ (current_gped[i] - initial_gped[i]) + 0.000000001
		);
		// ^^^^^ last small increment is to prevent any double-rounding fiasco's.
	}
}

void TFOOTMapProc::ProcessInitialPedestal() noexcept {
	if(N <= 0) ERROR("Starting to process, but local foot label <= 0 ?");

	const TFRSSortEvent* sortev = std::get<0>( this->in ).raw();
	
	auto [_0, Tlo, _1, _2, footN, data] = GetPtrs(sortev, N);
	if(footN == 0) return; 

	init_timing.assign(Tlo);
	
	/* A small barrier in case `dt_veto` field is given and consecutive timings fall below the dt veto. */
	/* `dt_veto` is already in nanoseconds, just like the timing.increment. */
	if(!std::isnan(dt_veto) and init_timing.initialized_ and init_timing.increment < dt_veto)
		return;

	int iraw;	
	for(int i=0; i < N_STRIPS; ++i) {
		iraw = (is_swapped == CableSwapped::YES) ? ((i + N_STRIPS/2) % N_STRIPS) : i;
		out.h2_raw_tmp->Fill(i, data[iraw]);
	}
}

//#define CALC_OFFSET_FROM_MEDIAN

void TFOOTMapProc::ProcessEventPedestal() noexcept {
	TFRSGo4Cont & input  = std::get<0>( this->in );
	out.Clean();

	const TFRSSortEvent* sortev = &input.inner();

	auto [_0, Tlo, _1, _2, footN, data] = GetPtrs(sortev, N);
	if(footN == 0) return;
	
	Scaler<32>& timing = out.inner().timing;
	timing.assign(Tlo);
	
	/* A small barrier in case `dt_veto` field is given and consecutive timings fall below the dt veto. */
	if(!std::isnan(dt_veto) and timing.initialized_ and timing.increment < dt_veto)
		return;
	
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
	
	for(int asic=0; asic < N_ASIC; ++asic) /* 0..=9 */ {
		const int i0 = asic * N_STRIPS_PER_ASIC;
		int i, iraw; 
		
		for(int strip=0; strip < N_STRIPS_PER_ASIC; ++strip) /* 0..=63 */ {
			i = i0 + strip; // <--- 'true strip' : 0..=640
			iraw = (is_swapped == CableSwapped::YES) ? ((i + N_STRIPS/2) % N_STRIPS) : i;
			
			out.h2_raw_tmp->Fill(i, data[iraw]);
			out.h2_raw    ->Fill(i, data[iraw]);

			/* Calculate the systematic shift per asic (group of 64 strips). */
			ped_offset[strip] = data[iraw] - current_gped[i];
		}

		/* Initial strip #0 is uncoupled. It contributes nothing. 
		 * Just remove it from the analysis. */
		ped_offset[0] = -DBL_MAX;

		std::sort(ped_offset.begin(), ped_offset.end());

		ped_off_med = mnd::median( ped_offset ); 

		ped_off_avg = std::accumulate (
			ped_offset.cbegin() + N_TRIM_FINE_PED_LO, 
			ped_offset.cend()   - N_TRIM_FINE_PED_HI, (double)0.0
		);
		ped_off_avg /= (N_STRIPS_PER_ASIC - N_TRIM_FINE_PED_LO - N_TRIM_FINE_PED_HI);

		out.h2_ped_off_med->Fill(asic, ped_off_med);
		out.h2_ped_off_avg->Fill(asic, ped_off_avg);
		out.h2_ped_off_diff->Fill(asic, ped_off_avg - ped_off_med);
	
		/* Repeat the loop, subtract the collective small offset (fine correction). */
		for(int strip=0; strip < N_STRIPS_PER_ASIC; ++strip) /* 0..=63 */ {
			i = i0 + strip;
			iraw = (is_swapped == CableSwapped::YES) ? ((i + N_STRIPS/2) % N_STRIPS) : i;

			adc_final = data[iraw] - current_gped[i];
			
			if(strip > 0) { /* The initial first strip in ASIC is uncoupled, just let it be. */
				adc_final -=
#ifdef CALC_OFFSET_FROM_MEDIAN
				ped_off_med
#else
				ped_off_avg
#endif
				; 
			} 
			else { /* Just for the initial strip that's uncoupled, wash away binning all the values into a single bin. */
				adc_final += rand() / (double)RAND_MAX ;
			}

			out.h2_corr->Fill(i, adc_final);
			out.inner().FOOTE[i] = adc_final;
		}
		out.inner().common_offset[asic] = 
#ifdef CALC_OFFSET_FROM_MEDIAN
			ped_off_med
#else
			ped_off_avg
#endif
		;
	}

	out.h1_entry_dt->Fill( out.inner().DeltaT().value_or(NAN) / 1000.0 /* in microseconds. */ );

	/* If we sample enough events, recalculate the global pedestal quickly from the batch. */
	if((++nsampled) == n_batch_pedestal.amount) {
		std::lock_guard<std::mutex> lock(g_foot_fit_mtx);
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
		parse_json_as_int_vec(parsed, *it);
	
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
		WARN("Ran over the data, but found 0 events with calibrated data? Setting all pedestals ridiculously high." 
			EMPH(FOOT: %d\n), N);
		for(int i=0; i < N_STRIPS; ++i) out.ped_s->at(i) = 10'000;
		return;
	}
	
	for(int i=0; i < N_STRIPS; ++i) {
		auto slice = std::unique_ptr<TH1D>( out.h2_corr->ProjectionY(mnd::msg("_py%d-%d", N, i), i+1, i+1) );
		slice->SetDirectory(nullptr);

		if(slice->GetEntries() < 100) {
			WARN("FOOT%d, strip %d, found no entries in y-projection of calibrated data? Skipping.\n", N, i);
			out.ped_s->at(i) = 10'000'000; /* Random number for high pedestal. */
			continue;
		}
		
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 10);
		double fitMax = slice->GetBinCenter(maxBin + 10);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);

		TF1* fitF = slice->GetFunction("gaus"); // Quiet null if no fit performed for empty data.

		if(fitF == nullptr) {
			WARN("FOOT%d, strip %d, couldn't fit Gauss for y-projection of calibrated data? Skipping. "
				"Debug: maxBin = %d, integral = %.1f\n", N, i, maxBin, slice->GetEntries());
			out.ped_s->at(i) = BAD_STRIP_CUTOFF_HI + 10;
		} else {
			out.ped_s->at(i) = fitF->GetParameter(2);
		}
	
		if(out.ped_s->at(i) > BAD_STRIP_CUTOFF_HI or 
			out.ped_s->at(i) < BAD_STRIP_CUTOFF_LO) {
			out.bad_strips->push_back(i);
		}
	}

	ParseStaticBadStrips();
}

/* ============= Extra aux JSON functions. ================= */

void TFOOTMapProc::parse_json_string(std::vector<int>& out, std::string s) {
	static const std::regex re_num(
		R"(^(0x|0b)?(\d+)$)"
	);
	static const std::regex re_range(
		R"(^(\d+)\.\.(=)?(\d+)$)"
	);
	static const std::regex re_seq(
		R"(^(\d+)n(\+\d+)?$)"
	);

	Trim(s);

	/* Strings can be passed either as:
	 * 1) raw numbers
	 * 2) range x1..x2 (x2 excluded) or x1..=x2 (x2 included)
	 * 3) sequence: an+b ('a', 'b' are the parameters, '+b' optional; 'n' fixed token). 
	 *    Meaning: strips: b, a+b, 2*a+b, etc. */
	std::smatch m;
	if(std::regex_match(s, m, re_num)) {
		if(m[1].matched) {
			if(m[1].str() == "0x") out.push_back(std::stoi(m[2].str(), nullptr, 16));
			else out.push_back(std::stoi(m[2].str(), nullptr, 2));
		}
		else out.push_back(std::stoi(m[2].str()));
	} else if(std::regex_match(s, m, re_range)) {
		int a = std::stoi(m[1].str());
		int b = std::stoi(m[3].str());
		if(a > b) ERROR("%d < %d found while parsing json: \'%s\'\n", a,b,s.c_str());
		for(int i=a; i<b; ++i) out.push_back(i);
		if(m[2].matched) out.push_back(b);
	} else if(std::regex_match(s, m, re_seq)) {
		int a = std::stoi(m[1].str());
		int b = m[2].matched ? std::stoi(m[2].str()) : 0;
		for(int i=b; i < 640; i+=a)
			out.push_back(i);
	} else {
		WARN("String \'%s\' doesn't match a: number, range or sequence regular expression.", s.c_str());
	}
}

/**
 * Json Custom range parser, accepting raw numbers or strings (or arrays of the same)
 * parsable as int or as a range:
 * 'a1..a2' left-inclusive, or 'a1..=a2' right inclusive.
 * */
void TFOOTMapProc::parse_json_as_int_vec(std::vector<int>& out, const json& j) {
	if(j.is_string()) {
		parse_json_string(out, j.get<std::string>());
	}
	else if(j.is_number()) {
		out.push_back(j.get<int>());	
	}
	else if(j.is_array()) {
		out.reserve(j.size());
		for(const auto& jsub : j) 
			parse_json_as_int_vec(out, jsub);
	}
	else ERROR("Passed a json object '%s' which isn't: string/array/number.\n", j.dump().c_str());
}

void TFOOTMapProc::append_flat_json(json& dst, const json& src) {
	if(!dst.empty() && !dst.is_object())
		ERROR("Destination object \'%s\' cannot store sources appended to it. Non-empty and not-an-object!", dst.dump().c_str());
	if(!src.is_object())
		ERROR("Source json: \'%s\' must be an object.", src.dump().c_str());

	for(const auto& [k, v] : src.items()) {
		auto it = dst.find(k);
		if(it == dst.end() || it->is_null()) {
			dst[k] = v;
			continue;
		}

		if(!it->is_array()) {
			*it = json::array({ *it }); // [old]
		}

		if(v.is_array()) {
			it->insert(it->end(), v.begin(), v.end()); // "k": [..., v[0], v[1], v[n-1] ]
		} else {
			it->push_back(v); // "k": [..., v]
		}
	}
}
