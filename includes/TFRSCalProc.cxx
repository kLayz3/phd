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
	for(auto& list : candidate_list)
		list.reserve(CANDIDATE_LIST_CAPACITY);
	for(auto& list : full_candidate_list)	
		list.reserve(CANDIDATE_LIST_CAPACITY * mnd::len(candidate_list));
}

void TFRSCalProc::ProcessEntry() noexcept {
	for(int i=0; i < N_VALID_SCI; ++i)
		this->ProcessSci(i);
	for(int i=0; i < N_VALID_TPC; ++i) {
		PreProcessTPC(i);
		for(int d : {0,1})
			ProcessDelayLine(i, d);		
		PostProcessTPC(i);
	}

	//this->ProcessS2Angle();
	//this->ProcessS4Angle();
}

void TFRSCalProc::PreProcessTPC(int _i_tpc) noexcept {
	RNTPCCal& out = (this->out).inner().tpc[_i_tpc];
	out.Clean();
}

/* `i` is the delay-line index: 0 or 1. */
void TFRSCalProc::ProcessDelayLine(int _i_tpc, int i) noexcept {
	constexpr auto INVALID = RNTPCMap::Measurement::TDC_INVALID;

	for(int a : {0,1} ) {
		candidate_list[a].clear();
		full_candidate_list[a].clear();
	}

	const auto& csum_lim    = TFRSCalCont::_tpc_param[_i_tpc].csum_lim;
	const auto& sci_ref_lim = TFRSCalCont::_tpc_param[_i_tpc].sci_ref_lim;

	const RNTPCMap& in = std::get<0>( this->in ).inner().tpc[_i_tpc];

	/* Input hits for this delay line. */
	const std::vector<RNTPCMap::Measurement>& hits = in.tdc[i];
	
	i32 csum = 0, csum_tmp;
	i32 ai = 0; 
	i32 tdc_a;
	
	/* Try to build up all the hits that satisfy the csum limits. */ 
	for(auto hit_l = hits.begin(); hit_l != hits.end() and hit_l->tdc_l != INVALID; ++hit_l) {
		const i32 dl = hit_l->tdc_l;

		for(auto hit_r = hits.begin(); hit_r != hits.end() and hit_r->tdc_r != INVALID; ++hit_r) {
			const i32 dr = hit_r->tdc_r;

			csum_tmp = dl + dr;

			for(int a : {0,1} ) {
				ai = 2*i + a; // csum limit array goes 0..4 not [2][2]
				
				for(auto hit_a = hits.begin(); hit_a != hits.end() and hit_a->tdc_a[a] != INVALID; ++hit_a) {
					tdc_a = hit_a->tdc_a[a];

					csum = csum_tmp - 2*tdc_a;	
					
					if(csum > csum_lim[ ai ][0] and csum < csum_lim[ ai ][1]) {
						candidate_list[a].emplace_back (
							tdc_a, dl, dr, INVALID	
						);
					}
				}
			}
		}
	}

	/* For this specific delay-line we collected all the hits that fit in a csum
	 * cut. Next do the y- selection which is much less selective than the csum one. 
	 * The testcase I kept failing was:
	 * DL(0): { [L: 40220, R: 40118, A: [33032, 33003]],
				[L: 45297, R: 46147, A: [38568, 38575]] }
	   DL(1): { [L: 40436, R: 40056, A: [33050, 33064]]
				[L: 45486, R: 46081, A: [38613, 38506]] }
		   S: { 21553, 28248 }
	 *
	 * In this case, since the csum selection is the strictest, it should correspond to a real
	 * particle. Since both 'S' values satisfy the y-check.
	 * Therefore, try to couple the anode TDC value to the lowest referent one, that fits in the cut.
	 * Having the map-level arrays sorted per-column comes in handy here!
	 */
	
	/* [1] Go over existing candidate hits, but sort them over the  `tdc_a` attribute.
	 * [2] There could be duplicate measurements here. Take the one with lower `tdc_l + tdc_r` score.  */
	for(auto& list : candidate_list) {
		if(list.size() < 2) continue;

		std::sort(list.begin(), list.end(),  // [1]
				[](const auto& lhs, const auto& rhs) { return lhs.a_tdc < rhs.a_tdc; });
		
		/* loop exit condition iterator is valid since vector at any point inside the loop must be sized >= 2. */
		for(auto it = list.begin(); it != list.end() - 1; ) { // [2]
			if( it->a_tdc == (it+1)->a_tdc ) {
				WARN("Found a duplicate in just one anode: TPC%d delay-line: %d\n", _i_tpc, i);
				auto next = it+1;
				if(it->dr_tdc + it->dl_tdc > next->dr_tdc + next->dl_tdc)
					it = list.erase(it); // Will point to `next`, but we don't assign directly due to cxxABI
				else
					list.erase(next), ++it; // Don't reassign, because `it` is still valid, it's before `next`.
			}
			else {
				++it;
			}
			/* list.end() iterator gets re-requested. */
		}
	}

	/* Try to see if any of the ref-sci values, if existent, fit in the y-cut. 
	 * The first one that does, we associate it with the hit and kick it out of the selection for
	 * consequent hits search. */
	
	auto it_end = in.tdc_ref.cend();
	
	for(int a : {0,1} ) { // This will just get unrolled twice, I can bet. :b
		auto& list = candidate_list[a];
		auto ref_it  = in.tdc_ref.cbegin(); // Start from begining, for a fresh anode.
		ai = 2*i + a;
		
		for(auto& hit : list) {
			for( ; ref_it != it_end; ++ref_it) {
				i32 diff = hit.a_tdc - *ref_it;
				if(diff > sci_ref_lim[ai][0] and diff < sci_ref_lim[ai][1]) {
					hit.ref_tdc = *ref_it++;
					break;
				}
			}
		}
	}
	
	/* On this delay line the corresponding hit list is complete.
	 * We can now turn it into actual position measurement: (x,y).
	 * Each of the two anode lists: candidate_list[0/1] are now also sorted according to
	 * ref tdc time.
	 * Combine their measurements. Alert if the same `ref_tdc` was found, but different dl/dr values.
	 * These must be consistent, otherwise we associate ghost particles. 
	 * If we find an anode candidate *without* matching `ref_tdc`, it gets thrown away, but maybe 
	 * similar candidate from the other anode managed to create a full (x,y) hit. */
	
	i32 current_ref_tdc = INVALID;
	auto it0 = candidate_list[0].cbegin();
	auto it1 = candidate_list[1].cbegin();
	
	/* Walk over anode(0); keep the corresponding iterator pointed to most relevant in anode(1) list. */
	for(; it0 < candidate_list[0].cend() and it0->ref_tdc != INVALID; ++it0) {
		current_ref_tdc = it0->ref_tdc;
		bool matched = false;

		/* In anode(1) list keep going until we encounter a hit with either
		 * [1] hit.ref_tdc == current_ref_tdc
		 *     => Combine the xy-measurements
		 * [2] hit.ref_tdc  > current_ref_tdc
		 *     => Save the iterator, maybe it matches on the next `it0`. */
		while( it1 < candidate_list[1].cend() and it1->ref_tdc != INVALID ) {
			if(it1->ref_tdc > current_ref_tdc) {
				/* Wait until anode(0) reference catches up. */
				break; 
			}

			else if(it1->ref_tdc < current_ref_tdc) {
				/* `it1` value must be unpaired, save it directly into the complete list. */
				full_candidate_list[i].emplace_back (
					std::array<i32,2>{ INVALID, it1->a_tdc },
					it1->dl_tdc,
					it1->dr_tdc,
					it1->ref_tdc
				);
			}

			else {
				/* The ref_tdc values are the same for both anodes. Delay line values must also match. */
				if(it1->dl_tdc != it0->dl_tdc or it1->dr_tdc != it0->dr_tdc) {
					/* This measurement is just broken then somehow. Kick it out. */
					WARN("TPC%d, delay-line: %d. Found two VALID hits without matching delay-line measurement?\n"
						"        A     L     R     S\n"
						"0: (%4d %4d %4d %4d)\n1: (%4d %4d %4d %4d)\n",
						_i_tpc, i, 
						it0->a_tdc, it0->dl_tdc, it0->dr_tdc, it0->ref_tdc,
						it1->a_tdc, it1->dl_tdc, it1->dr_tdc, it1->ref_tdc);
					++it1;
				}
				else {	
					matched = true;	
				}
			}
		}

		if(matched) {
			full_candidate_list[i].emplace_back (
				std::array<i32,2>{ it0->a_tdc, it1->a_tdc },
				it0->dl_tdc,
				it0->dr_tdc,
				current_ref_tdc
			);
		} else {
			full_candidate_list[i].emplace_back (
				std::array<i32,2>{ it0->a_tdc, INVALID },
				it0->dl_tdc,
				it0->dr_tdc,
				current_ref_tdc
			);
		}
	}

	/* Walk over possible remnants in anode(1); anode(0) has been iterated through. */
	for(; it1 < candidate_list[1].cend() and it1->ref_tdc != INVALID; ++it1)
		full_candidate_list[i].emplace_back (
			std::array<i32,2>{ INVALID, it1->a_tdc },
			it1->dl_tdc,
			it1->dr_tdc,
			it1->ref_tdc
		);


}

void TFRSCalProc::PostProcessTPC(int _i_tpc) noexcept {
	constexpr auto INVALID = RNTPCMap::Measurement::TDC_INVALID;
	RNTPCCal& out = (this->out).inner().tpc[_i_tpc];

	const auto& [bx, ax, by, ay, _, __, ___] = TFRSCalCont::_tpc_param[_i_tpc]; // read from .rodata
	for(int d : {0,1} ) {
		auto& l = full_candidate_list[d];
		std::sort(l.begin(), l.end());

		std::vector<RNTPCCal::Measurement>& out_list = out.hits[d];
		out_list.reserve(l.size());
		
		double x,y;
		int a0 = 2*d;
		for(const auto& hit : l) {
			int mask = 0;

			/* x-pos: */ 
			x = ax[d] * ( hit.dl_tdc - hit.dr_tdc)  + bx[d];
			/* y-pos: */ 
			y = [&]() noexcept -> double {
				double tmp = 0;
				if(hit.a_tdc[0] != INVALID) {
					tmp += ay[a0] * ( hit.a_tdc[0] - hit.ref_tdc) + by[a0];
					mask |= 1;
				}
				if(hit.a_tdc[1] != INVALID) {
					tmp += ay[a0+1] * ( hit.a_tdc[1] - hit.ref_tdc) + by[a0+1];
					mask |= 2;
				}
				if(mask == 3) tmp /= 2;
				return tmp;
			}();

			out_list.emplace_back(x,y, hit.ref_tdc, mask);
		}
	}	
}

#if 0
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

#endif
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
