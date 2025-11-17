#include "core/AuxFunctions.hh"

#include "TFRSCalProc.h"
#include "TFRSMapCont.h"

constexpr auto INVALID = RNTPCMap::Measurement::TDC_INVALID;

TFRSCalProc::TFRSCalProc(TFRSCalCont& out, const TFRSMapCont& in) : TFRSCalProc::Base(out, in) {
	tpc_hits.reserve(HIT_LIST_CAPACITY);
	for(auto& anode_list : candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
	for(auto& anode_list : initial_candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
	
	for(int i=0; i<N_VALID_TPC; ++i) {
		if(! out.tpc_param) ERROR("TPC parameter field not constructed?\n");
		tpc_param[i] = (*out.tpc_param)[i];
	}
}

void TFRSCalProc::ProcessEntry() noexcept {
	ncalled++;
	for(int i=0; i < N_VALID_TPC; ++i)
		this->ProcessTPC(i);				
}

void TFRSCalProc::ProcessTPC(int _i_tpc) noexcept {	
	tpc_hits.clear();
	for(auto& c : candidate_list) c.clear();
	for(auto& c : initial_candidate_list) c.clear();
	
	const TFRSMapCont& input = std::get<0>( this->in );
	const RNTPCMap& in = input.inner().tpc[_i_tpc];
	auto const& hits = in.tdc; 

	RNTPCCal& out = this->out.inner().tpc[_i_tpc];

	/* For a general multihit combination, try to see which falls into Control Sum limits.
	 * Ideally a hit lights up a delay line (LR), corresponding 1-2 anodes and the referent scintillator. */
	
	const auto& [bx, ax, by, ay, csum_lim, sci_ref_lim] = tpc_param[_i_tpc]; 
	
	/* Go over all recorded hits from map. Different channels can have different multihit value.
	 * They're padded with `INVALID` value to indicate it's missing from a channel.
	 * It's asserted that in the list of hits, for any channel, datum values will always come before
	 * invalid fill. */
	 
	 /* In general, the following algorithms are asymptotically bad - n*m*k complexity, 
	 * but these multihits are very small on the order of ~2-4. 
	 * Proper binary search, instead of brute force, would be faster starting from sizes 10x10x10 or so. 
	 * Because each of the 9 sequences of RNTPCMap::Measurement 
	 * {tdc_l[0]}, {tdc_l[1]}, {tdc_r[0]}, are sorted until (optional) INVALID range. */

	/* First selection - only select anode(i) and sci_ref values which can have their Y-reconstructed. */
	for(const auto& h0 : hits) {
		if(h0.tdc_ref == INVALID) break;
		
		for(const auto& h1 : hits) {
			/* Try to see if any of the anode values, if existent, fit in the y-cut. */
			for(int a=0; a<4; ++a) {
				if(h1.tdc_a[a] != INVALID) {
					auto diff = h1.tdc_a[a] - h0.tdc_ref;
					if(diff > sci_ref_lim[a][0] and diff < sci_ref_lim[a][1])
						initial_candidate_list[a].emplace_back(h1.tdc_a[a], INVALID, INVALID, h0.tdc_ref);
				}
			}
		}
	}
	
	/* Second selection - go over above collected anode hits, 
	 * and select those that satisfy their corresponding control sum check,
	 * with some pair (l,d) of delay-line measurement. */

	int csum1, csum2, ref; int _cval;
	for(int a=0; a<4; ++a) {
		ref = a >> 1;
		const auto& [lim_x_lo, lim_x_hi] = csum_lim[a];
		
		for(auto& candidate : initial_candidate_list[a]) {
			_cval = (candidate.a_tdc << 1);
			for(const auto& h1 : hits) { // Delay-line left potential pairing partners
				if(h1.tdc_l[ref] == INVALID) break;
				csum1 = h1.tdc_l[ref] - _cval; 

				for(const auto& h2 : hits) { // Delay-line right potential pairing partners
					if(h2.tdc_r[ref] == INVALID) break;

					csum2 = h2.tdc_r[ref] + csum1;
					if(csum2 > lim_x_hi) break; /* Following entries will have even higher `tdc_r` */

					/* Found a valid candidate! */
					else if(csum2 > lim_x_lo) { 
						candidate.dl_tdc = h1.tdc_l[ref];	
						candidate.dr_tdc = h1.tdc_r[ref];
						if( IsUniqueTPCMeasurement(candidate_list[a], candidate) )
							candidate_list[a].push_back(candidate);
					}
				}
			}
		}
	}

	/* Write the values to the output. */

	/* [1] Find how many elements we need to allocate. */
	size_t _n_measurements = 0;
	for(int a=0; a<4; ++a)
		_n_measurements = std::max (
			candidate_list[a].size(),
			_n_measurements
		);

	out.hits.resize(_n_measurements);

	for(int a=0; a<4; ++a) {
		auto ref = a >> 1;
		for(int n=0; n < (int)candidate_list[a].size(); ++n) {
			TPCHitCandidate& valid_hit = candidate_list[a][n];
			out.hits.at(n).at(a) = {
				/* x-pos: */ ax[ref] * ( valid_hit.dl_tdc - valid_hit.dr_tdc)  + bx[ref],
				/* y-pos: */ ay[ a ] * ( valid_hit.a_tdc  - valid_hit.ref_tdc) + by[ a ]
			};
		}
	}
}

bool TFRSCalProc::IsUniqueTPCMeasurement(const TPCHitCandidateList& list, const TPCHitCandidate& candidate) noexcept {
	/* Candidate must be completely unique to qualify a good measurement.
	 * Do a small linear search to see if it can be added. */
	for(const auto& e : list) 
		if(e.a_tdc   == candidate.a_tdc  ||
		   e.dl_tdc  == candidate.dl_tdc ||
		   e.dr_tdc  == candidate.dr_tdc ||
		   e.ref_tdc == candidate.ref_tdc
		) return false;
	return true;
}
