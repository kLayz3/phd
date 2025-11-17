#include "TFOOTCalProc.h"
#include "TFOOTMapCont.h"

#include "TGraph.h"

#include <algorithm>
#include <cmath>

template<> bool TFOOTCalProc::_IsAddibleToCluster<TFOOTCalProc::kPOS>(const int );
template<> bool TFOOTCalProc::_IsAddibleToCluster<TFOOTCalProc::kNEG>(const int );

bool operator==(TFOOTCalProc::ClusterFullType& full, TFOOTCalProc::ClusterType t) {
	using E = TFOOTCalProc::ClusterType;
	auto& l = full.first;
	auto& r = full.second;
	
	switch(t) {
		case(E::kGOOD):
			return (l == E::kGOOD) and (r == E::kGOOD);
		case(E::kFRAGMENTED):
			return (
				(r != E::kMERGED) and
				(l == E::kFRAGMENTED or r == E::kFRAGMENTED)
			);
		case(E::kWAVEY):
			return (
				(r != E::kMERGED and l != E::kMERGED) and
				(r == E::kWAVEY or l == E::kWAVEY)
			);
		case(E::kMERGED):
			return (r == E::kMERGED or l == E::kMERGED);
		default:
			return (l == E::kUNKNOWN or r == E::kUNKNOWN);
	}
}

bool operator!=(TFOOTCalProc::ClusterFullType& full, TFOOTCalProc::ClusterType t) {
	return !(full == t) ; 
}

TFOOTCalProc::ClusterType TFOOTCalProc::GetClusterType() {
	using E = TFOOTCalProc::ClusterType;
	if(_ct == E::kGOOD)       return E::kGOOD;
	if(_ct == E::kMERGED)     return E::kMERGED;
	if(_ct == E::kFRAGMENTED) return E::kFRAGMENTED;
	if(_ct == E::kWAVEY)      return E::kWAVEY;
	
	return E::kUNKNOWN;
}

/* Input and output should be properly assigned a priori. */
TFOOTCalProc::TFOOTCalProc(TFOOTCalCont& out, TFOOTMapCont& in, double x_seed, double x_neigh) : 
	TFOOTCalProc::Base(out, in), 
	X_CENTRE_THR(x_seed), X_NEIGHB_THR(x_neigh)
{
	TFOOTMapCont& input = std::get<0>(this->in);
	if(!input.ped_s || input.ped_s->size() != N_STRIPS)
		ERROR("Input sigma array not initialized? State (null|entries): \'%s\' . Did you call TFOOTMapCont::Setup() ?\n",
				input.ped_s ? "null" : Form("%zu", input.ped_s->size()));

	if(out.FOOT_N < 0 || out.POS < 0) 
		ERROR("Output object (TFOOTCalCont): \'%s\' uninitialized. Did you call ::Setup?", out.GetName());

	const auto& sigma = *input.ped_s;
	for(int i=0; i<N_STRIPS; ++i) {
		c_thr[i] = sigma[i] * X_CENTRE_THR; 
		n_thr[i] = sigma[i] * X_NEIGHB_THR; 
		//fprintf(stderr, "%.2f ", sigma[i]);
	} //fprintf(stderr, "\n");

	/* Bad/dead strips just label them with large thresholds.. */
	const std::vector<int>& bad_strips = *input.bad_strips;
	for(auto i : bad_strips) {
		c_thr[i] = BAD_STRIP_FAKE_THRESHOLD; 
		n_thr[i] = BAD_STRIP_FAKE_THRESHOLD; 
	}

	for(int i=0; i<N_STRIPS; ++i) {
		if(n_thr[i] < 0.01)
			ERROR("IDK bro: FOOT%d, strip = %d, n_thr = %.2f. Fix it.\n", input.FOOT_N, i, n_thr[i]);
		if(c_thr[i] < 0.01)
			ERROR("IDK bro: FOOT%d, strip = %d, c_thr = %.2f. Fix it.\n", input.FOOT_N, i, c_thr[i]);
	}

	/* Set the _e pointer. */
	_e = &input.inner().FOOTE[0];

#define _ALTER_TITLE(x) \
	x->SetTitle(Form("%s : S=%.1f N=%1.f", x->GetTitle(), X_CENTRE_THR, X_NEIGHB_THR))
	_ALTER_TITLE(out.h1_cl_type);
	_ALTER_TITLE(out.h1_raw_mult);
	_ALTER_TITLE(out.h1_mult);
	_ALTER_TITLE(out.h1_X);
	_ALTER_TITLE(out.h1_dE);
	_ALTER_TITLE(out.h1_sn_ratio);
}

void TFOOTCalProc::ProcessEntry() noexcept {
	out.Clean();
	_e = &std::get<0>(in).inner().FOOTE[0]; /* Don't know if rebinding is really necessary... */

	if( std::isnan(_e[0]) ) return; /* Missing data; in previous step marked it NAN by default. */
	
	/* Copy	the data over, since there's a write to `e`. */
	memcpy(e, _e, sizeof(e)); 
	int i = 0;
	/* Try to find a central 'seed' strip. */
	for(; i < N_STRIPS; ++i) {
		if(e[i] > c_thr[i])
			MakeACluster(i);
	}
}

void TFOOTCalProc::MakeACluster(int& c0 /* Starting index. Passes C-threshold check. */ ) {
	_ClearClust();

	/* From the loop caller, it's assumed that previous two 
	 * strips aren't part of this specific cluster. A.K.A, they're below the N-threshold. 
	 * And that the current strip `i` is above the C-threshold. */
	
	ClusterType& _ct_neg = _ct.first;  /* Negative side. */
	ClusterType& _ct_pos = _ct.second; /* Positive side. */

	const int i_init = c0;
	
	int i = i_init;
	_AddHit(i);
	
	/* The moment `e[i] > n_thr[i]` fails, a lookahead to i+1, or a lookbehind to i-1 occurs.
	 * If the lookahead/behind finds a hit, e[i] is averaged out between e[i-1] and e[i+1].
	 * If it fails, means both e[i] and:
	 *     > e[i+1], in case of right traversing, 
	 *     > e[i-1], in case of left traversing,
	 * fail the N-threshold cutoff and it constitutes a 
	 * cluster end, for this respective side's traverse.  */

	while((++i) < N_STRIPS and (e[i] > n_thr[i] || _IsAddibleToCluster<kPOS>(i)))  {
		_AddHit(i);
	}
	/* From here, `i` points to end of cluster. */
	c0 = i;

	/* In case this positive (right) side connected a fragment, restore the 'old' energy. */
	if(_ct_pos == ClusterType::kFRAGMENTED) {
		auto [_i, _e] = _frag;
		this->e[_i] = _e;
	}

	i = i_init;
	while((--i) >= 0 and (e[i] > n_thr[i] || _IsAddibleToCluster<kNEG>(i))) {
		_AddHit(i);
	}

	/* In case this negative (left) side connected a fragment, restore the 'old' energy. */
	if(_ct_neg == ClusterType::kFRAGMENTED) {
		auto [_i, _e] = _frag;
		this->e[_i] = _e;
	}
	
	/* Handle cluster size == 1.
	 * Most of the time it's noise or some other weird effect.
	 * However, possible that indeed there's a just-below-threshold effect for small
	 * central hits. */
	if(_cl_cnt == 1) {
		/* Bad charge sharing can only be for MIP's, or if the 
		 * central strip somehow crossed a 'bad strip'. */
		if(_buf[0].e > 3 * c_thr[i]) return;
	}

	/* Check if cluster is 'merged'. Meaning that in the cluster,
	 * there exist two strips above C-threshold with 2 or more strips in-between 
	 * which are below C-threshold. */
	/* Sort by their position. */
	std::sort(_buf, _buf + _cl_cnt, [](const auto& l, const auto& r) { return l.x < r.x; });

	TClustHit* ref_hit = nullptr;
	for(u32 i=0; i < _cl_cnt; ++i) {
		auto [strip, adc_val] = _buf[i];
		if( adc_val > c_thr[strip] ) {
			if(ref_hit && (strip - ref_hit->x > 2)) {
				_ct_pos = ClusterType::kMERGED;
				break;
			}
			ref_hit = &_buf[i];
		}
	}

	/* Handle the wavyness in the cluster. If it is fragmented/merged, then it's wavy by default.
	 * Meaning: first sequence of collected strips energy must be monotonically increasing,
	 * while second sequence must have monotonically decreasing energy values, as to yield a proper
	 * hit structure, if the cluster is to be marked `kGOOD`.  */
	
	if(_ct != ClusterType::kFRAGMENTED and _ct != ClusterType::kMERGED) {
		double prev_e = 0;
		enum { rising, falling } curve = rising;

		for(u32 hit = 0; hit < _cl_cnt; ++hit) {
			double curr_e = _buf[hit].e;
			if(curve == rising and curr_e < prev_e) {
				curve = falling;
			}
			else if(curve == falling and curr_e > prev_e) {
				_ct_pos = ClusterType::kWAVEY;
				break;
			}
			prev_e = curr_e;
		}
	}
	
	/* Integrate energy, and make a weighted average for position. */
	double cl_e = 0, cl_wx = 0;
	for(u32 hit = 0; hit < _cl_cnt; ++hit) {
		cl_e  += _buf[hit].e;
		cl_wx += _buf[hit].e * _buf[hit].x;
	}
	
	cl_wx /= cl_e;

	auto [cl_max_i, cl_max_e] = *std::max_element(_buf, _buf + _cl_cnt, 
		[](const auto& l, const auto& r) { 
			return l.e < r.e; 
		});
	double cl_m = cl_e / cl_max_e;

	ClusterType ct = this->GetClusterType();
	out.inner().AddCluster(cl_wx, cl_e, cl_m, ct);
	out.h1_raw_mult->Fill(_cl_cnt);
	out.h1_mult->Fill(cl_m);
	out.h1_dE->Fill(cl_e);
	out.h1_X->Fill(cl_wx);
	out.h1_cl_type->Fill(ct);

	if(_cl_cnt >= MASSIVE_CLUSTER_CUTOFF && out.inner()._fBadE.size() == 0) {
		out.inner()._fBadE.assign(e, e + N_STRIPS); 
	}

	switch(_cl_cnt) {
		case 1:
			out.h1_dE_m1->Fill(cl_e);
			break;
		case 2:
			out.h1_dE_m2->Fill(cl_e);
			break;
		case 3:
			out.h1_dE_m3->Fill(cl_e);
			break;
		default:
			break;
	}
	/* Since the cluster hit array is sorted by the strip index, and only consecutive 
	 * sequence of hits makes a cluster - means this below is valid: */
	cl_max_i -= _buf[0].x;
	/* assert(cl_max_i >= 0 && cl_max_i < _cl_cnt); */

	if(cl_max_i > 0 && cl_max_i+1 < (int)_cl_cnt) {
		double e_left  = _buf[cl_max_i - 1].e;
		double e_right = _buf[cl_max_i + 1].e;
		double sn = 1 / log( cl_max_e * cl_max_e / (e_left*e_right) );
		out.h1_sn_ratio->Fill(sn);
	}

	if(cl_max_e > 30 and _cl_cnt == 1 and 
		out.inner()._fHeClSize1.size() == 0 and 
		out.inner()._fBadE.size() == 0) 
	{
		out.inner()._fHeClSize1.assign(e, e + N_STRIPS);
	}
}

template<>
bool TFOOTCalProc::_IsAddibleToCluster<TFOOTCalProc::kPOS>(const int i /* 1, ... , 639 */) {
	/* Lookahead by one strip.
	 * It's known from caller's logic that `i-1` belongs to the cluster and
	 * that `i` is below the N-threshold. */

	ClusterType& _ct_pos = _ct.first; /* Positive side. */

	/* Don't allow multiple fragmentations of cluster. A true zig-zag pattern. */
	if(_ct_pos == ClusterType::kFRAGMENTED) return false;

	int prev = i-1 /* <-- Is in cluster. */ , next = i+1;
	if(next < N_STRIPS and e[next] > n_thr[next]) {
		_frag = {i, e[i]};
		_ct_pos = ClusterType::kFRAGMENTED;
		
		e[i] = ( e[prev] + e[next] ) / 2;

		// Calling function will restore this hit, here just 'reassign' the energy.	
		return true; 
	}
	return false;
}

template<>
bool TFOOTCalProc::_IsAddibleToCluster<TFOOTCalProc::kNEG>(const int i /* 0, ... , 638 */) {
	/* Lookbehind by one strip.
	 * It's known from caller's logic that `i+1` belongs to the cluster and
	 * that `i` is below the N-threshold. */

	ClusterType& _ct_neg = _ct.first; /* Negative side. */

	/* Don't allow multiple fragmentations of cluster. A true zig-zag pattern. */
	if(_ct_neg == ClusterType::kFRAGMENTED) return false;
	
	int next = i+1 /* <-- Is in cluster. */ , prev = i-1;
	if(prev >= 0 and e[prev] > n_thr[prev]) {
		_frag = {i, e[i]};
		_ct_neg = ClusterType::kFRAGMENTED;
		
		e[i] = ( e[prev] + e[next] ) / 2;
		
		// Calling function will restore this hit, here just 'reassign' the energy.
		return true; 
	}
	return false;
}
