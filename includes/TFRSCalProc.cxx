#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"
#include "Eigen/Dense"

/* This ostream API is really hacked. Close your eyes. */
#ifdef TFRSCALPROC_VERBOSE_
template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& );

std::ostream& operator<<(std::ostream& os, const RNTPCMap::Measurement& m) {
	os << "{" << KBH_MAG << "A: " << KNRM << m.tdc_a;
	os  << KBH_MAG << " L: " << KNRM << m.tdc_l 
		<< KBH_MAG << " R: " << KNRM << m.tdc_r  << "}";
	return os;
}

#include "helper_fwd.h"

std::ostream& operator<<(std::ostream& os, const RNTPCMap& m) {
	os << KBH_RED << "Full TPC dump: " << std::endl;
	os << KBH_GRN << "DL0: " << KNRM << m.tdc[0] << '\n' 
		<< KBH_GRN << "DL1: " << KRNM << m.tdc[1] << '\n';
	os << KBH_BLU << "S:   " << KRNM << BOLD << m.tdc_ref << KNRM;
	return os;
}
#endif
/* Now code becomes human again... */

/* Compare if the first argument lands inside the area of two delimiters. */
template<typename T, typename U>
bool is_inside(T x, const std::array<U,2>& lim) {
	return (lim[0] <= x) and (x <= lim[1]);
}

/* Return a number in interval: [0,1> */ 
static double uniform() noexcept {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1ULL); 
}

TFRSCalProc::TFRSCalProc(TFRSCalCont& out, const TFRSMapCont& in) : TFRSCalProc::Base(out, in) {
	for(auto& l2 : candidate_list)
		for(auto& list : l2) list.reserve(CANDIDATE_LIST_CAPACITY);
	for(auto& list : full_candidate_list)	
		list.reserve(CANDIDATE_LIST_CAPACITY * mnd::len(candidate_list));
}

void TFRSCalProc::ProcessEntry() noexcept {
#ifdef TFRSCALPROC_VERBOSE_
	++called_;
#endif
	for(int i=0; i < N_VALID_SCI; ++i)
		this->ProcessSci(i);
	for(int i=0; i < N_VALID_TPC; ++i)
		this->ProcessTPC(i);

	//this->ProcessS2Angle();
	//this->ProcessS4Angle();
}

void TFRSCalProc::ProcessTPC(int _i_tpc) noexcept {
	PreProcessTPC(_i_tpc);
	const RNTPCMap& in = std::get<0>( this->in ).inner().tpc[_i_tpc];

	/* If the reference signal is empty, no way to reconstruct y. Just skip this TPC. */
	if(in.tdc_ref.size() == 0) return;

	for(int d : {0,1} /* dl index */) 
		ProcessCSum(_i_tpc, d);

#ifdef TFRSCALPROC_SINGLEHIT
	if(in.tdc_ref.size() != 1) return;
	i32 ref_tdc = in.tdc_ref[0];
	
	ProcessSingleHit(_i_tpc, ref_tdc);
#else
	for(int d : {0,1} /* dl index */) 
		ResolveTPCDelayLineConflicts(_i_tpc, d);
	
	FilterRemainingConflicts();

	/* At this point, each CSum (tuple: {a_tdc, l_tdc, r_tdc}) was found and marked as
	 * valid or invalid. Next construct the y-position. */
	ProcessTPC_YCut(_i_tpc);
#endif
	PostProcessTPC(_i_tpc);
}

void TFRSCalProc::PreProcessTPC(int _i_tpc) noexcept {
	RNTPCCal& out = (this->out).inner().tpc[_i_tpc];
	out.Clean();

	/* Reset the local buffers to process this TPC. */
	for(auto& dl_list : candidate_list)
		for(auto& anode_list : dl_list) anode_list.clear();

	for(auto& dl_list : conflicts)
		for(auto& anode_list : dl_list) anode_list.clear();

	for(int d: {0,1})
		full_candidate_list[d].clear();
}

/* `i` is the delay-line index: 0 or 1. */
void TFRSCalProc::ProcessCSum(int _i_tpc, int i) noexcept {
	constexpr i32 INVALID = RNTPCMap::Measurement::TDC_INVALID;
	
	const RNTPCMap& in = std::get<0>( this->in ).inner().tpc[_i_tpc];

	const auto& csum_lim = TFRSCalCont::_tpc_param[_i_tpc].csum_lim; // read from .rodata
	auto& dl_list = candidate_list[i];

	/* Input hits for this delay line. */
	const std::vector<RNTPCMap::Measurement>& hits = in.tdc[i];
	
	i32 csum, csum_tmp0, csum_tmp1;
	i32 ai;
	i32 tdc_a, dl, dr;
	
	/* Try to build up **all** the hits that satisfy the csum limits. */ 
	for(int a : {0,1} ) {
		auto& list = dl_list[a]; 
		ai = 2*i + a; // csum limit array goes 0..4 not [2][2]
		const std::array<double,2>& csum_limits = csum_lim[ai];
	
		/* Loop over all anode[a] hits. */
		for(auto hit_a = hits.begin(); hit_a != hits.end() and hit_a->tdc_a[a] != INVALID; ++hit_a) {
			tdc_a = hit_a->tdc_a[a];
			csum_tmp0 = -2 * tdc_a;

			/* Independent loop over all dl_l hits. */
			for(auto hit_l = hits.begin(); hit_l != hits.end() and hit_l->tdc_l != INVALID; ++hit_l) {
				dl = hit_l->tdc_l;
				csum_tmp1 = csum_tmp0 + dl;

				/* Independent loop over all dl_r hits. */
				for(auto hit_r = hits.begin(); hit_r != hits.end() and hit_r->tdc_r != INVALID; ++hit_r) {
					dr = hit_r->tdc_r;
					csum = csum_tmp1 + dr;	

					if( is_inside(csum, csum_limits) ) {
						list.emplace_back (
							tdc_a, dl, dr, INVALID	
						);
					}
				}
			}
		}
	} 

	/* The testcase I kept failing was:
	 * DL(0): { [L: 40220, R: 40118, A: [33032, 33003]],
				[L: 45297, R: 46147, A: [38568, 38575]] }
	   DL(1): { [L: 40436, R: 40056, A: [33050, 33064]]
				[L: 45486, R: 46081, A: [38613, 38506]] }
		   S: { 21553, 28248 }
	 *
	 * In this case, since the csum selection is the strictest, it should correspond to a real
	 * particle. Since both 'S' values satisfy the y-check.
	 * Therefore, try to couple the anode TDC value to the lowest referent one, that fits in the cut.
	 * Having the map-level arrays sorted per-column comes in handy here! */
	
	/* Principle is: we can reconstruct completely **ONLY** as many hits as there are sci- referent hits. */

#ifndef TFRSCALPROC_SINGLEHIT
	/* There could be conflicting measurements already here. Just mark their iterator. */
	for(int a : {0,1} ) {
		auto& list = dl_list[a]; 
		if(list.size() < 2) continue;
			
		ai = 2*i + a;

		/* loop exit condition iterator is valid since vector at any point inside the loop must be sized >= 2. */
		for(auto it0 = list.begin(); it0 != list.end(); ++it0) {
			for(auto it1 = it0+1; it1 != list.end(); ++it1) {
			if( it0->a_tdc  == it1->a_tdc || 
				it0->dl_tdc == it1->dl_tdc ||
				it0->dr_tdc == it1->dr_tdc ) {
			
				it0->status = TFRSCalProc::TPCHitCandidate::Status::kIN_CONFLICT; 
				it1->status = TFRSCalProc::TPCHitCandidate::Status::kIN_CONFLICT;

				conflicts[i][a].emplace_back(&*it0, &*it1);

#ifdef TFRSCALPROC_VERBOSE_
					WARN(EMPH([%lu])": Found a CSUM duplicate in just one anode: TPC%d delay-line: %d, anode %d.\n"
						"           A     L     R\n"
						"Hit(0): (%4d %4d %4d)\nHit(1): (%4d %4d %4d)\n",
						called_, _i_tpc, i, a, 
						it0->a_tdc, it0->dl_tdc, it0->dr_tdc,
						it1->a_tdc, it1->dl_tdc, it1->dr_tdc);
					std::cerr << in << "\nCsumLim: " << csum_lim[ai] << "\n\n";
					++errors_;
#endif
				}
			}
		}
	}
#endif

}

#ifndef TFRSCALPROC_SINGLEHIT
void TFRSCalProc::ResolveTPCDelayLineConflicts(int _i_tpc, int i /* dl id: 0,1 */) noexcept { 
	std::array<TPCConflicts, 2>& conflicts_in_dl = conflicts[i];

	const auto& anode_diff_lim = TFRSCalCont::_tpc_param[_i_tpc].anode_diff_lim[i];
	const auto& dl_left_diff_lim = TFRSCalCont::_tpc_param[_i_tpc].dl_left_diff_lim;
	const auto& dl_right_diff_lim = TFRSCalCont::_tpc_param[_i_tpc].dl_right_diff_lim;
	
	for(int a : {0,1} /* anode id: 0,1 */) {
		TPCConflicts& conflict = conflicts_in_dl[a];

		for(auto& cf : conflict) {
			/* If either party got dropped in previous conflict, mark the other one
			 * as resolved and move on. */
			if(cf.first->status == TPCHitCandidate::Status::kDROPPED) {
				cf.second->status = TPCHitCandidate::Status::kFINE;
				continue;
			} else if(cf.second->status == TPCHitCandidate::Status::kDROPPED) {
				cf.first->status = TPCHitCandidate::Status::kFINE;
				continue;
			}

			/* Possibly both `IN_CONFLICT` or one got resolved 'fine' from previous resolutions.
			 * Either way, do the test again. */
			bool sgn_value = static_cast<bool>(i);
			if( TryResolveViaDiff<MeasurementType::kLEFT>(cf, dl_left_diff_lim, candidate_list[i^1][0], sgn_value) ) continue;
			if( TryResolveViaDiff<MeasurementType::kRIGHT>(cf, dl_right_diff_lim, candidate_list[i^1][0], sgn_value) ) continue;
			if( TryResolveViaDiff<MeasurementType::kLEFT>(cf, dl_left_diff_lim, candidate_list[i^1][1], sgn_value) ) continue;
			if( TryResolveViaDiff<MeasurementType::kRIGHT>(cf, dl_right_diff_lim, candidate_list[i^1][1], sgn_value) ) continue;
			
			/* Last resort, very unlikely to succeed. */
			if( TryResolveViaDiff<MeasurementType::kANODE>(cf, anode_diff_lim, candidate_list[i][a^1], static_cast<bool>(a) ) ) continue;
		}
	}
}

void TFRSCalProc::FilterRemainingConflicts() noexcept { 
	for(auto& dl_line_lists : candidate_list) {
		for(auto& list : dl_line_lists) {
			mnd::Erase(list, [](const auto& c) { return c.status != TFRSCalProc::TPCHitCandidate::Status::kFINE; } );
		}
	}
}

/* Try and resolve conflict `c` by hints from the list of other anode/dl. 
 * 's' parameter is to switch the signum of the difference, based on if the argument is number 0 or 1.
 * In parameters, it is always [0]-[1] . So if the signum is '1' the sign swap needs to happen. */
template<TFRSCalProc::MeasurementType m>
bool TFRSCalProc::TryResolveViaDiff(TPCConflict& c, const std::array<double,2>& limits, const TPCHitCandidateList& other, bool s) {
	std::pair<double, double> diff;
	std::pair<u16,u16> wins {0,0};

	for(const auto& o : other) {
		u8 result = 0;
		
		if constexpr(m == MeasurementType::kLEFT) {
			diff.first  = c.first ->dl_tdc - o.dl_tdc;
			diff.second = c.second->dl_tdc - o.dl_tdc;
		}
		else if constexpr(m == MeasurementType::kRIGHT) {
			diff.first  = c.first ->dr_tdc - o.dr_tdc;
			diff.second = c.second->dr_tdc - o.dr_tdc;
		}
		else if constexpr(m == MeasurementType::kANODE) {
			diff.first  = c.first ->a_tdc - o.a_tdc;
			diff.second = c.second->a_tdc - o.a_tdc;
		}
		if(s) { diff.first = -diff.first; diff.second = -diff.second; };

		if( is_inside(diff.first, limits) ) {
			result |= 1;
		} else if( is_inside(diff.second, limits) ) {
			result |= 2;
		}
		
		switch(result) {
			case 1:
				++wins.first;
				break; 
			case 2: 
				++wins.second;
				break;
			default: break; /* A tie. */
		}
	}
	if(wins.first > wins.second) {
		c.first ->status = TPCHitCandidate::Status::kFINE;
		c.second->status = TPCHitCandidate::Status::kDROPPED;
		return true;
	} else if(wins.second > wins.first) {
		c.first ->status = TPCHitCandidate::Status::kDROPPED;
		c.second->status = TPCHitCandidate::Status::kFINE;
		return true;
	}

	return false;
}

void TFRSCalProc::ProcessTPC_YCut(int _i_tpc) noexcept { 
	const RNTPCMap& in = std::get<0>( this->in ).inner().tpc[_i_tpc];
	const u32 N = (u32)in.tdc_ref.size();
	
	const auto& ylim = TFRSCalCont::_tpc_param[_i_tpc].sci_ref_lim;
	/* Per anode, we can expect at most as many valid hits as there are references.
	 * This is usually 1-2 in the scintillator. */

	/* TPC hit candidates (per anode) are naturally sorted according to `a_tdc` field. */
	for(int d: {0,1}) {
		auto& list = candidate_list[d];
		/* Now we try to bundle up the 2 anode measurements 
		 * `list[0]` and `list[1]` the same. */
			
	}
}
#else
void TFRSCalProc::ProcessSingleHit(int _i_tpc, int ref_tdc) noexcept {
	constexpr i32 INVALID = RNTPCMap::Measurement::TDC_INVALID;

	for(int d: {0,1}) {
		auto& outlist = full_candidate_list[d];
		
		for(int a: {0,1}) {
			auto& list = candidate_list[d][a];
			if(list.size() != 1) continue;
			
			if(outlist.empty()) {
				outlist.emplace_back (
					INVALID, INVALID,
					list[0].dl_tdc, list[0].dr_tdc,
					ref_tdc
				);
				outlist[0].a_tdc[a] = list[0].a_tdc;
			}
			else {
				/* Check if dl-dr measurements match. */
				if( outlist[0].dl_tdc == list[0].dl_tdc &&
					outlist[0].dr_tdc == list[0].dr_tdc )
				{
					outlist[0].a_tdc[a] = list[0].a_tdc;
				} 
				else { /* Somehow this anode's hit doesn't have matching dl/dr. Kick out the whole thing. */
#ifdef TFRSCALPROC_VERBOSE_
					static u64 called_1 = 0;
					if(_i_tpc == 0 and d == 0) 
						WARN(EMPH2([%lu])": Found a CSUM valid in two anode with mismatched measurements?: TPC%d delay-line: %d.\n"
							"           A     L     R\n"
							"Ano(%d): (%4d %4d %4d)\nAno(%d): (%4d %4d %4d)\n",
							++called_1, _i_tpc, d, 
							a, list[0].a_tdc, list[0].dl_tdc, list[0].dr_tdc,
							a^1, outlist[0].a_tdc[a^1], outlist[0].dl_tdc, outlist[0].dr_tdc);

#endif
					outlist.pop_back();
				}
			}
		}
	}
}
#endif

void TFRSCalProc::PostProcessTPC(int _i_tpc) noexcept {
	constexpr i32 INVALID = RNTPCMap::Measurement::TDC_INVALID;
	RNTPCCal& out = (this->out).inner().tpc[_i_tpc];

	[[ maybe_unused ]]
	const auto& [bx, ax, by, ay, _1, _2, _3, _4, _5, _6] = TFRSCalCont::_tpc_param[_i_tpc]; // read from .rodata
	
	for(int d : {0,1} ) {
		auto& l = full_candidate_list[d];
		//std::sort(l.begin(), l.end());
#ifdef TFRSCALPROC_SINGLEHIT
		assert(l.size() < 2 && mnd::msg("What, TPC%d, has %lu hits, and should be 1?", l.size()));
#endif
		std::vector<RNTPCCal::Measurement>& out_list = out.hits[d];
		out_list.reserve(l.size());
	
		int a0 = 2*d;
		for(const auto& hit : l) {
			/* hit.a_tdc[0/1] is nullable, handle it explicitly here. */
			double diff_tdc_0 = (hit.a_tdc[0] != INVALID) ? static_cast<double>( hit.a_tdc[0] - hit.ref_tdc ) : NAN; 
			double diff_tdc_1 = (hit.a_tdc[1] != INVALID) ? static_cast<double>( hit.a_tdc[1] - hit.ref_tdc ) : NAN; 

			out_list.emplace_back (
				/* x : */ ax[d] * ( hit.dl_tdc - hit.dr_tdc)  + bx[d],
				/* y0: */ ay[a0+0] * diff_tdc_0 + by[a0+0],
				/* y1: */ ay[a0+1] * diff_tdc_1 + by[a0+1],
				/* s : */ hit.ref_tdc,
				/* tr: */ 1
			);

			const RNTPCCal::Measurement& m = out_list.back();
			this->out.h2_tpc_xy[_i_tpc][d] -> Fill( m.X(), m.Y() );
			this->out.h1_tpc_mask[_i_tpc][d] -> Fill( m.AnodeMask() );
		}
	}
}

/* Preliminary angle analysis at S2. */
void TFRSCalProc::ProcessS2Angle() noexcept {
#if 0 
	constexpr int N = 3; /* Take 3 TPC's into the fit. */

	static const auto& tpc_param = this->out._tpc_param;
	static const std::array<double,4> z = {
		tpc_param[0].z0,
		tpc_param[1].z0,
		tpc_param[2].z0,
		tpc_param[3].z0,
	};

	constexpr double zT = 3355 - 440/2.0; 
	using QR = Eigen::ColPivHouseholderQR<Eigen::Matrix<double, N, 2>>;
	static QR qr = {};

	[[maybe_unused]]
	static const Eigen::Matrix<double, N, 2> A = []{
		Eigen::Matrix<double, N, 2> tmp;
		for(int i=0; i<N; ++i) {
			tmp(i,0)  = z[i];
			tmp(i,1) = 1.0;
		}
		qr.compute(tmp);
		return tmp;
	}();
	
	const auto& tpc = this->out.inner().tpc;
	if(tpc[0].hits[0].size() == 1 and tpc[0].hits[1].size() == 1 and
		tpc[1].hits[0].size() == 1 and tpc[1].hits[1].size() == 1 and 
		tpc[2].hits[0].size() == 1 and tpc[2].hits[1].size() == 1)
	{
		std::array<double, N> x = {0}, y = {0};
		for(int i=0; i<N; ++i) {
			for(int d : {0,1}) {
				x[i] += tpc[i].hits[d][0].x;
				y[i] += tpc[i].hits[d][0].y;
			}
			x[i] /= 2;
			y[i] /= 2;
		}
		Eigen::Map<const Eigen::Matrix<double, N, 1>> xv(x.data());
		Eigen::Map<const Eigen::Matrix<double, N, 1>> yv(y.data());

		Eigen::Vector2d coeffs_x = qr.solve(xv);
		Eigen::Vector2d coeffs_y = qr.solve(yv);
		double ax = coeffs_x(1), bx = coeffs_x(0);
		double ay = coeffs_y(1), by = coeffs_y(0);
		
		out.h2_ab_s2_before_target->Fill(1000*ax, 1000*bx);
		out.h2_xy_s2_at_target->Fill(ax*zT + bx, ay*zT + by); 
	}
	
	if(tpc[3].hits[0].size() == 1 and tpc[3].hits[1].size() == 1) {
		double x = (tpc[3].hits[0][0].x + tpc[3].hits[1][0].x) / 2;
		double y = (tpc[3].hits[0][0].y + tpc[3].hits[1][0].y) / 2;
		out.h2_xy_s2_after_target->Fill(x,y);
		out.h2_ab_s2_after_target->Fill (
			1000 * x / (z[3] - zT),
			1000 * y / (z[3] - zT)
		);
	}
#endif
}

/* Preliminary angle analysis at S2. */
void TFRSCalProc::ProcessS4Angle() noexcept {}

void TFRSCalProc::ProcessSci(int _i_sci) noexcept {
	constexpr i32 INVALID = RNSciMap::Measurement::TDC_INVALID;
	
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
