#include "TFOOTPedestalProc.h"
#include "TF1.h"
#include "libs.hh"
#include "dbg.hh"
#include <cassert>
#include <numeric>
#include "AuxFunctions.hh"

static_assert(TFOOTPedestalProc::N_STRIPS == TFOOTPedestalProc::N_ASIC * TFOOTPedestalProc::N_STRIPS_PER_ASIC, "Failed build: nstrip != nasic*nstrip_per_asic!\n");

void TFOOTPedestalProc::ProcessEntry() noexcept {
	switch(process_type) {
		case kGPED: 
			ProcessGlobalPedestal();
			break;
		case kEPED:
			ProcessEventPedestal();
			break;
	}
}

void TFOOTPedestalProc::ProcessGlobalPedestal() noexcept {
	if(*data._FOOT == 0) return;
	FOR(i, N_STRIPS) {
		data.h2_raw->Fill(i, data._FOOTE[i]);
	}
}

/* This is NOT thread safe to run in parallel since `slice->Fit` inherently talks to gROOT, gPad.
 * Even if TH1D histograms are detached. */
void TFOOTPedestalProc::CalcGlobalPedestal() {
	if(data.h2_raw->GetEntries() == 0) {
		dbg("Ran over the TTree, but found 0 events with data?", data.FOOT_N, "Setting all pedestals to 0.");
		return;
	} 
	FOR(i, N_STRIPS) {
		TH1D* slice = data.h2_raw->ProjectionY("", i+1, i+1);
		slice->SetDirectory(nullptr);
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 20);
		double fitMax = slice->GetBinCenter(maxBin + 20);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		data.gped->at(i)   = fitF->GetParameter(1);
		data.gped_s->at(i) = fitF->GetParameter(2);
		data.gr_s0->SetPoint(i, i, data.gped_s->at(i));
		
		delete slice;
	}
}

//#define CALC_OFFSET_FROM_MEDIAN

void TFOOTPedestalProc::ProcessEventPedestal() noexcept {
	this->Clear();
	if(*data._FOOT == 0) return;
	
	/* Subtract the global pedestal. */
	/* Algorithm (per groups of 64-strips) is the following: j=0,1,2, ... 63
	 * (1) Calculate the offsets of each: _FOOTE[j] - gped[j] and store in a 64-length array.
	 * (2) Sort the array.
	 * (3) Trim the lowest and highest N_TRIM number of entries
	 * (4) From what is remaining, (re)calculate the mean, yielding the offset `o[asic]` of each asic
	 * (5) Offset the calculation yielding the final calibrated ADC value:
	 * --> adc[i] = _FOOTE[i] - gped[i] - o[i / N_STRIPS_PER_ASIC] */

	/* Should be between 5..10 */
#define N_TRIM_FINE_PED 6

	double ped_off_med, ped_off_avg, adc_intermediate, adc_final; 
	FOR(asic, N_ASIC) /* 0..=9 */ {
		const int i0 = asic * N_STRIPS_PER_ASIC;
		int i; 
		
		/* Calculate the systematic shift per asic (group of 64 strips). */
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip;
			ped_offset[strip] = data._FOOTE[i] - data.gped->at(i);
		}
		std::sort(ped_offset.begin(), ped_offset.end());

		ped_off_med = median(ped_offset); 
		ped_off_avg = std::accumulate(ped_offset.cbegin() + N_TRIM_FINE_PED, ped_offset.cend() - N_TRIM_FINE_PED, (double)0.0);
		ped_off_avg /= (N_STRIPS_PER_ASIC - 2*N_TRIM_FINE_PED);

		data.h2_ped_off_med->Fill(asic, ped_off_med);
		data.h2_ped_off_avg->Fill(asic, ped_off_avg);
		data.h2_ped_off_diff->Fill(asic, ped_off_avg - ped_off_med);
	
		FOR(strip, N_STRIPS_PER_ASIC /* 0..=63 */) {
			i = i0 + strip;
			adc_intermediate = data._FOOTE[i] - data.gped->at(i);
			data.h2_mid->Fill(i, adc_intermediate);
			adc_final = adc_intermediate -
#ifdef CALC_OFFSET_FROM_MEDIAN
				ped_off_med
#else
				ped_off_avg
#endif
			;
			data.h2_corr->Fill(i, adc_final); 
			data.FOOTE[i] = adc_final;
		}
	}
}

void TFOOTPedestalProc::CalcFinalPedestal() {
	if(data.h2_corr->GetEntries() == 0) {
		dbg("Ran over the TTree, but found 0 events with calibrated data?", data.FOOT_N);
		return;
	} 
	FOR(i, N_STRIPS) {
		TH1D* slice = data.h2_corr->ProjectionY("", i+1, i+1);
		slice->SetDirectory(nullptr);
		int maxBin = slice->GetMaximumBin();
		double fitMin = slice->GetBinCenter(maxBin - 10);
		double fitMax = slice->GetBinCenter(maxBin + 10);

		slice->Fit("gaus", "Q", "", fitMin, fitMax);
		TF1* fitF = slice->GetFunction("gaus");
		data.gped_sf->at(i) = fitF->GetParameter(2);
		data.gr_s1->SetPoint(i, i, data.gped_sf->at(i));
		
		delete slice;
	}
}
