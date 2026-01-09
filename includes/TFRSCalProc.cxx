#include "Eigen/Core"
#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"
#include "Eigen/Dense"

/* Return a number in interval: [0,1> */ 
static double uniform() noexcept {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1ULL); 
}

TFRSCalProc::TFRSCalProc(TFRSCalCont& out, const TFRSMapCont& in) : TFRSCalProc::Base(out, in) {
	for(auto& anode_list : candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
	for(auto& anode_list : initial_candidate_list)
		anode_list.reserve(CANDIDATE_LIST_CAPACITY);
	
	full_candidate_list.reserve(CANDIDATE_LIST_CAPACITY * mnd::len(candidate_list));
}

void TFRSCalProc::ProcessEntry() noexcept {
	for(int i=0; i < N_VALID_SCI; ++i)
		this->ProcessSci(i);
	for(int i=0; i < N_VALID_TPC; ++i)
		this->ProcessTPC(i);

	this->ProcessS2Angle();
	this->ProcessS4Angle();
}

void TFRSCalProc::ProcessTPC(int _i_tpc) noexcept {	
	constexpr auto INVALID = RNTPCMap::Measurement::TDC_INVALID;

	for(auto& c : candidate_list) c.clear();
	for(auto& c : initial_candidate_list) c.clear();
	full_candidate_list.clear();
	
	const TFRSMapCont& input = std::get<0>( this->in );
	const RNTPCMap& in = input.inner().tpc[_i_tpc];
	auto const& hits = in.tdc; 

	RNTPCCal& out = this->out.inner().tpc[_i_tpc];

	/* For a general multihit combination, try to see which falls into Control Sum limits.
	 * Ideally a hit will light up a delay line (LR), corresponding 1 or both anodes, and the referent scintillator. */
	
	[[maybe_unused]]
	const auto& [bx, ax, by, ay, csum_lim, sci_ref_lim, _] = TFRSCalCont::_tpc_param[_i_tpc]; 
	
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
							h0.tdc_ref   /* Scintillator reference. */
						);
					}
				}
			} // loop over anodes
		} // loop over all hits; only anode parts interesting 
	} // loop over all hits; only SCI part interesting
	
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
				} // loop over potential candidates
			} // loop over initial candidates
		}
	} // loop over anodes indices 0,1,2,3

	/* For each anode (a = 0,1,2,3) the corresponding hit list is complete.
	 * But it's possible it isn't completely sorted. So just sort it now. 
	 * Sorting is done based on referent (Sci) TDC value. Check bool operator<(...) function in the header. */
	/* We extend the hitcandidate type to also attach the anode index with it. */

	/* [1] Create this extended hit list. */ 
	/* [2] Find how many elements we need to allocate,
	 * for each complete anode measurement(s). */
	for(int a=0; a<4; ++a) {
		const auto& list = candidate_list[a];	
		for(const auto& hit : list)
			full_candidate_list.emplace_back(a, hit);
	}
	/* [3] Sort the list w.r.t. sci reference TDC. */
	std::sort(full_candidate_list.begin(), full_candidate_list.end());
	
	/* [4] Calculate how many output elements we need to allocate. */
	size_t _n_measurements = 0;
	int curr_ref_tdc = -1; 
	for(const auto& hit_extended : full_candidate_list) {
		if(hit_extended.hit.ref_tdc != curr_ref_tdc) {
			++_n_measurements;
			curr_ref_tdc = hit_extended.hit.ref_tdc;
		}
	}
	//static u64 _n = 0;
	//WARN("Ev: [%lu] _n_measurements = %zu\n", _n, _n_measurements);

	/* [5] Resize output container. */
	out.hits.assign(_n_measurements, { /* Default ctor == nan all fields */ });

	/* [6] Sort out the hits into the output container. Do the linear calibration. */
	curr_ref_tdc = -1;
	int n = -1; 
	for(const auto& ehit : full_candidate_list) {
		const auto a   = ehit.index; /* Anode index: 0,1,2,3. */	
		const auto ref = a >> 1;     /* Corresponding delay line index: 0,1. */
		const TPCHitCandidate& valid_hit = ehit.hit;
		
		if(valid_hit.ref_tdc != curr_ref_tdc) {
			++n;
			curr_ref_tdc = valid_hit.ref_tdc;
		}
		out.hits.at(n).at(a) = {
			/* x-pos: */ ax[ref] * ( valid_hit.dl_tdc - valid_hit.dr_tdc)  + bx[ref],
			/* y-pos: */ ay[ a ] * ( valid_hit.a_tdc  - valid_hit.ref_tdc) + by[ a ],
		};
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

/* This gets checked only per anode. */
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

/* Preliminary angle analysis at S2. */
void TFRSCalProc::ProcessS2Angle() noexcept {
	constexpr int n = 3;

	static const auto& tpc_param = this->out._tpc_param;
	static const std::array<double, n> z = {
		tpc_param[0].z0,
		tpc_param[1].z0,
		tpc_param[2].z0
	};
	
	static const Eigen::Matrix<double, n, 2> A = [] {
		Eigen::Matrix<double,n,2> tmp;
		for(int i=0; i<n; ++i) {
			tmp(i,0)  = z[i];
			tmp(i,1) = 1.0;
		}
		return tmp;
	}();

	const auto& tpc = this->out.inner().tpc;
	if(tpc[0].hits.size() == 1 and
		tpc[1].hits.size() == 1 and
		tpc[2].hits.size() == 1 )
	{
		
		double x[n] = {0}, y[n] = {0};
		for(int i=0; i<n; ++i) {
			for(int a=0; a<4; ++a) {
				x[i] += tpc[i].hits[0][a].x;
				y[i] += tpc[i].hits[0][a].y;
			}
			x[i] /= 4;
			y[i] /= 4;
		}
		Eigen::VectorXd X(n);
		Eigen::VectorXd Y(n);

		for(int i=0; i<n; ++i) {
			X(i)    = x[i];
			Y(i)    = y[i];
		}
		Eigen::Vector2d coeffs_x = A.colPivHouseholderQr().solve(X);
		Eigen::Vector2d coeffs_y = A.colPivHouseholderQr().solve(Y);
		double a = coeffs_x(0);
		double b = coeffs_y(0);
		
		out.h2_ab_s2_before_target->Fill(1000*a, 1000*b);
	}

}

/* Preliminary angle analysis at S2. */
void TFRSCalProc::ProcessS4Angle() noexcept {
	
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

	[[maybe_unused]]
	const auto& [bx, ax, lim, _] = TFRSCalCont::_sci_param.at(_i_sci);

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
