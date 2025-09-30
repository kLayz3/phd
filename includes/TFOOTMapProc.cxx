#include "TFOOTMapProc.h"
#include "TF1.h"
#include "libs.hh"
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <numeric>
#include <unordered_set>
#include "AuxFunctions.hh"

static_assert(TFOOTMapProc::N_STRIPS == TFOOTMapProc::N_ASIC * TFOOTMapProc::N_STRIPS_PER_ASIC, "Failed build: nstrip != nasic*nstrip_per_asic!\n");

using nlohmann::json;
json TFOOTMapProc::_bad_strips{};

void TFOOTMapProc::ProcessEntry() noexcept {
	switch(process_type) {
		case kGPED: 
			ProcessGlobalPedestal();
			break;
		case kEPED:
			ProcessEventPedestal();
			break;
	}
}

void TFOOTMapProc::ProcessGlobalPedestal() noexcept {
	if(*data->_FOOT == 0) return;
	int iraw;
	FOR(i, N_STRIPS) {
		iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;
		data->h2_raw->Fill(i, data->_FOOTE[iraw]);
	}
}

/* This is NOT thread safe to run in parallel since `slice->Fit` inherently talks to gROOT, gPad.
 * Even if TH1D histograms are detached. */
void TFOOTMapProc::CalcGlobalPedestal() {
	if(data->h2_raw->GetEntries() == 0) {
		WARN("Ran over the TTree, but found 0 events with data?" EMPH(FOOT: %d) ", Setting all pedestals to 0.", data->FOOT_N);
		return;
	} 
	FOR(i, N_STRIPS) {
		TH1D* slice = data->h2_raw->ProjectionY("", i+1, i+1);
		slice->SetDirectory(nullptr);
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 20);
		double fitMax = slice->GetBinCenter(maxBin + 20);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		data->gped->at(i)   = fitF->GetParameter(1);
		data->gped_s->at(i) = fitF->GetParameter(2);
		data->gr_s0->SetPoint(i, i, data->gped_s->at(i));
		
		delete slice;
	}
}

//#define CALC_OFFSET_FROM_MEDIAN

void TFOOTMapProc::ProcessEventPedestal() noexcept {
	data->Clean();
	if(*data->_FOOT == 0) return;
	
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

	double ped_off_med, ped_off_avg, adc_intermediate, adc_final; 
	FOR(asic, N_ASIC) /* 0..=9 */ {
		const int i0 = asic * N_STRIPS_PER_ASIC;
		int i; 
		
		/* Calculate the systematic shift per asic (group of 64 strips). */
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip;
			int iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;
			ped_offset[strip] = data->_FOOTE[iraw] - data->gped->at(i);
		}
		/* Initial strip #0 is uncoupled. It contributes nothing. 
		 * Just remove it from the analysis. */
		ped_offset[0] = -DBL_MAX;

		std::sort(ped_offset.begin(), ped_offset.end());

		ped_off_med = median(ped_offset); 
		ped_off_avg = std::accumulate(ped_offset.cbegin() + N_TRIM_FINE_PED_LO, ped_offset.cend() - N_TRIM_FINE_PED_HI, (double)0.0);
		ped_off_avg /= (N_STRIPS_PER_ASIC - N_TRIM_FINE_PED_LO - N_TRIM_FINE_PED_HI);

		data->h2_ped_off_med->Fill(asic, ped_off_med);
		data->h2_ped_off_avg->Fill(asic, ped_off_avg);
		data->h2_ped_off_diff->Fill(asic, ped_off_avg - ped_off_med);
	
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip;
			int iraw = is_swapped ? ((i + N_STRIPS/2) % N_STRIPS) : i;
			adc_intermediate = data->_FOOTE[iraw] - data->gped->at(i);
			data->h2_mid->Fill(i, adc_intermediate);
			adc_final = adc_intermediate;
			if(strip > 0) {
				adc_final -=
#ifdef CALC_OFFSET_FROM_MEDIAN
				ped_off_med
#else
				ped_off_avg
#endif
				; /* The first strip in ASIC is uncoupled, just let it be. */
			} 
			else { /* For the initial strip that's uncoupled, wash away binning all the values into a single bin, . */
				adc_final += rand() / (double)RAND_MAX ;
			}
			data->h2_corr->Fill(i, adc_final); 
			data->inner().FOOTE[i] = adc_final;
		}
	}
}

void TFOOTMapProc::CalcFinalPedestal() {
	if(data->h2_corr->GetEntries() == 0) {
		WARN("Ran over the TTree, but found 0 events with calibrated data?" EMPH(FOOT: %d\n), data->FOOT_N);
		return;
	}
	FOR(i, N_STRIPS) {
		TH1D* slice = data->h2_corr->ProjectionY("", i+1, i+1);
		slice->SetDirectory(nullptr);
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 10);
		double fitMax = slice->GetBinCenter(maxBin + 10);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		data->gped_sf->at(i) = fitF->GetParameter(2);
		data->gr_s1->SetPoint(i, i, data->gped_sf->at(i));
	
		if(fitF->GetParameter(2) > BAD_STRIP_CUTOFF_HI or 
			fitF->GetParameter(2) < BAD_STRIP_CUTOFF_LO) {
			data->bad_strips->push_back(i);
		}
		delete slice;
	}

	this->ParseStaticBadStrips();
}

int TFOOTMapProc::ParseStaticBadStrips() {
	if(!data || !data->bad_strips)
		ERROR("output bad strips container uninitialized. Did you call TFOOTMapCont::Init(..) beforehand?");

	const char* key = Form("FOOT%d", data->FOOT_N);
	std::vector<int> parsed{};
	if(auto it = _bad_strips.find(key); it != _bad_strips.end())
		::parse_json_as_int_vec(parsed, *it);
	
	auto& v = *data->bad_strips;
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

