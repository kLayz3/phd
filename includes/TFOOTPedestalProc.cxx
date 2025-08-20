#include "TFOOTPedestalProc.h"
#include "TF1.h"
#include "libs.hh"
#include "dbg.hh"
#include <cassert>
#include <numeric>
#include "AuxFunctions.hh"

static_assert(TFOOTPedestalProc::N_STRIPS == TFOOTPedestalProc::N_ASIC * TFOOTPedestalProc::N_STRIPS_PER_ASIC, "Failed build: nstrip != nasic*nstrip_per_asic!\n");

void TFOOTPedestalProc::ProcessEntry() noexcept {
	if(do_global_pedestal) ProcessGlobalPedestal();
	else ProcessEventPedestal();
}

void TFOOTPedestalProc::ProcessGlobalPedestal() noexcept {
	if(data._FOOT == 0) return;
	FOR(i, N_STRIPS)
		data.h2_raw->Fill(i, data._FOOTE[i]);
}

void TFOOTPedestalProc::CalcGlobalPedestal() {
	if(data.h2_raw->GetEntries() == 0) {
		dbg("Ran over the TTree, but found 0 events with data?", data.FOOT_N, "Setting all pedestals to 0.");
		return;
	}
	FOR(i, N_STRIPS) {
		TH1D* slice = data.h2_raw->ProjectionY(Form("_FOOT_PROJ%d", i+1), i+1);
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 6);
		double fitMax = slice->GetBinCenter(maxBin + 6);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		data.gped->at(i) = fitF->GetParameter(1);
		data.gped_s->at(i) = fitF->GetParameter(2);
		data.h2_s0->Fill(i, data.gped_s->at(i));
	}
#if 0
	FOR(asic, N_ASIC) {
		gped_sum_per_asic[asic] = std::accumulate(
				data.gped->cbegin() + asic * N_STRIPS_PER_ASIC,
				data.gped->cbegin() + (asic + 1) * N_STRIPS_PER_ASIC,
				0.0
			);
	}
#endif
}

//#define CALC_OFFSET_FROM_MEDIAN

void TFOOTPedestalProc::ProcessEventPedestal() noexcept {
	this->Clear();
	if(data._FOOT == 0) return;
	
	data.has_data = true;

	/* Subtract the global pedestal. */
	/* Algorithm (per groups of 64-strips) is the following:
	 * (1) Calculate the offsets of each _FOOTE[i] - gped[i] and store in an array.
	 * (2) Sort the array. 
	 * (3) Trim the lowest and highest N_TRIM number of entries
	 * (4) from what's remaining, (re)calculate the mean, yielding the offset `o[asic]` of each asic
	 * (5) Offset the calculation yielding the final calibrated ADC value:
	 * --> adc[i] = _FOOTE[i] - gped[i] - o[i / N_STRIPS_PER_ASIC] */

	/* Should be between 5..10 */
#define N_TRIM_FINE_PED 6

	double ped_off_med, ped_off_avg, adc_final; 
	FOR(asic, N_ASIC) {
		const int i0 = asic * N_STRIPS_PER_ASIC;
		int i; 
		
		/* Calculate the systematic shift per asic. */
		FOR(strip, N_STRIPS_PER_ASIC /* 0..64 */) {
			i = i0 + strip;
			ped_offset[strip] = data._FOOTE[i] - data.gped->at(i);
		}
		std::sort(ped_offset.begin(), ped_offset.end());
		ped_off_med = median(ped_offset); 
		ped_off_avg = std::accumulate(ped_offset.cbegin() + N_TRIM_FINE_PED, ped_offset.cend() - N_TRIM_FINE_PED, 0);
		ped_off_avg /= (N_STRIPS_PER_ASIC - 2*N_TRIM_FINE_PED);

		data.h2_ped_off_med->Fill(asic, ped_off_med);
		data.h2_ped_off_avg->Fill(asic, ped_off_avg);
		data.h2_ped_off_diff->Fill(asic, ped_off_avg - ped_off_med);
	
		FOR(strip, N_STRIPS_PER_ASIC /* 0..64 */) {
			i = i0 + strip;
			adc_final = data._FOOTE[i] - data.gped->at(i) -
#ifdef CALC_OFFSET_FROM_MEDIAN
				ped_off_med
#else
				ped_off_avg
#endif
			;
			data.FOOTE[i] = adc_final;
		}
	}
}
