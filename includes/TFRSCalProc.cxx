#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"

/* Return a number in interval: [0,1> */ 
double uniform() noexcept {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1ULL); 
}

TFRSCalProc::TFRSCalProc(TFRSCalCont& out, const TFRSMapCont& in) : TFRSCalProc::Base(out, in) {
	for(auto& anode_list : candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
	for(auto& anode_list : initial_candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
}

void TFRSCalProc::ProcessEntry() noexcept {
	for(int i=0; i < N_VALID_SCI; ++i)
		this->ProcessSci(i);
	for(int i=0; i < N_VALID_TPC; ++i)
		this->ProcessTPC(i);
}

void TFRSCalProc::ProcessTPC(int _i_tpc) noexcept {	
	constexpr auto INVALID = RNTPCMap::Measurement::TDC_INVALID;

	for(auto& c : candidate_list) c.clear();
	for(auto& c : initial_candidate_list) c.clear();
	
	const TFRSMapCont& input = std::get<0>( this->in );
	const RNTPCMap& in = input.inner().tpc[_i_tpc];
	auto const& hits = in.tdc; 

	RNTPCCal& out = this->out.inner().tpc[_i_tpc];

	/* For a general multihit combination, try to see which falls into Control Sum limits.
	 * Ideally a hit will light up a delay line (LR), corresponding 1 or both anodes, and the referent scintillator. */
	
	const auto& [bx, ax, by, ay, csum_lim, sci_ref_lim] = TFRSCalCont::_tpc_param[_i_tpc]; 
	
	/* Go over all recorded hits from map. Different channels can have different multihit value.
	 * They're padded with `INVALID` value to indicate it's missing from a channel.
	 * It's asserted that in the list of hits, for any channel, datum values will always come before
	 * invalid fill. */
	 
	/* In general, the following algorithms are asymptotically bad - n*m*k complexity, 
	 * but these multihits are very small on the order of ~2-4. 
	 * Proper binary search, instead of brute force, would be faster starting from sizes 10x10x10 or so. 
	 * Because each of the 9 sequences of RNTPCMap::Measurement 
	 * {tdc_l[0]}, {tdc_l[1]}, {tdc_r[0], ...}, are sorted until (optional) INVALID range. */

	/* First selection - go over all `sci_ref` values and only select anode(i) which can have their Y-reconstructed. */
	for(const auto& h0 : hits) {
		if(h0.tdc_ref == INVALID) break;
		
		for(const auto& h1 : hits) {
			/* Try to see if any of the anode values, if existent, fit in the y-cut. */
			for(int a=0; a<4; ++a) {
				if(h1.tdc_a[a] != INVALID) {
					auto diff = h1.tdc_a[a] - h0.tdc_ref;
					if(diff > sci_ref_lim[a][0] and diff < sci_ref_lim[a][1]) {
						initial_candidate_list[a].emplace_back(
							h1.tdc_a[a], /* Anode TDC value. */
							INVALID,     /* Delay-left  TDC (to-be-found). */
							INVALID,     /* Delay-right TDC (to-be-found). */ 
							h0.tdc_ref /* Scintillator reference. */
						);
					}
				}
			}
		}
	}
	
	/* Second selection - go over all four anodes (a), and the corresponding
	 * candidate hit list: (a, _, _, tdc_ref), 
	 * and select pairs (l,d) of corresponding delay-line measurements, 
	 * that satisfy their corresponding control sum check. */

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

					else if(csum2 > lim_x_lo) {  
						/* Found a valid candidate! */
						candidate.dl_tdc = h1.tdc_l[ref];	
						candidate.dr_tdc = h1.tdc_r[ref];
						
						if( IsUniqueTPCMeasurement(candidate_list[a], candidate) )
							candidate_list[a].push_back(candidate);
					}
				}
			}
		}
	}

	/* For each anode (a = 0,1,2,3) the corresponding hit list is complete.
	 * But it's possible it isn't completely sorted. So just sort it now. 
	 * Sorting is done based on referent (Sci) TDC value. */
	for(int a=0; a<4; ++a) {
		auto& list = candidate_list[a];	
		std::sort(list.begin(), list.end());
	}

	/* Write the values to the output. */

	/* [1] Find how many elements we need to allocate,
	 * for each complete anode measurement(s). */
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
				/* y-pos: */ ay[ a ] * ( valid_hit.a_tdc  - valid_hit.ref_tdc) + by[ a ],
				/* ref  : */ valid_hit.ref_tdc
			};
		}
	}
	
	/* Make a preliminary plot of (x,y), only if _n_measurements is 1. */
	if((_i_tpc == 2 or _i_tpc == 3) and _n_measurements == 1) {
		const auto& hit = out.hits[0];
		double x0=0, y0=0;
		int n=0;
		for(const auto& anode : hit) 
			if(! std::isnan(anode.x)) {
				x0 += anode.x; y0 += anode.y;
				++n;
			}
		x0 /= n;
		y0 /= n;
		if(_i_tpc == 2) this->out.h2_xy_s2_before_target->Fill(x0, y0);
		else if(_i_tpc == 3) this->out.h2_xy_s2_after_target->Fill(x0, y0);
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

void TFRSCalProc::ProcessSci(int _i_sci) noexcept {
	constexpr auto INVALID = RNSciMap::Measurement::TDC_INVALID;
	
	RNSciMap const& in = std::get<0>(this->in).inner().sci[_i_sci];
	RNSciCal&      out = (this->out).inner().sci[_i_sci];
	out.Clean();
	
	out.E = static_cast<double>(in.qdc[0] * in.qdc[1] + 0.0000001);
	out.E = sqrt(out.E);
	/* ^^^ Last addition is t make sqrt() stable for (0,0) combination. */

	const auto& hits = in.tdc;

	const auto& [bx, ax, lim] = TFRSCalCont::_sci_param.at(_i_sci);

	double d_l, d_r;

	/* Try to associate left and right hits. */
	for(int i=0; i < (int)hits.size(); ++i) {
		if(hits[i].tdc_l == INVALID) break;

		for(int j=i; (j<(int)hits.size()) && (hits[j].tdc_r != INVALID); ++j) {
			auto diff_lr = hits[i].tdc_l - hits[j].tdc_r;
			if(diff_lr > lim[0] and diff_lr < lim[1]) {
				/* Draw twice to account for TDC uncertainty. */
				d_l = uniform(); d_r = uniform();
				auto avg_t = (hits[i].tdc_l + hits[j].tdc_r + d_l + d_r) / 2;
				double x = ax * (diff_lr + d_l - d_r) + bx;

				double t = SCIParam::channel_to_ns * avg_t;
				out.hits.emplace_back(x,t);
				
				break;
			}
		}
	}

	if(_i_sci == 0) {
		for(const auto& [x,t] : out.hits)
			this->out.h1_x_sc21_before_target->Fill(x);
	}
	if(_i_sci == 1) {
		for(const auto& [x,t] : out.hits)
			this->out.h1_x_sc22_after_target->Fill(x);
	}
}
