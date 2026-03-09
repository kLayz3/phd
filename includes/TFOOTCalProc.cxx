#include "TFOOTCalProc.h"
#include "TFOOTCalCont.h"
#include "TFOOTMapCont.h"
#include "TH2I.h"

#include <algorithm>
#include <cmath>

#include "json_struct_def.hh"
#include "FastGauss.hxx"

// This will make all the fits be analytical from 3 points around the maximum
// intensity strip.
//#define FIT_ALL_WITH_3_POINTS

/* Fit function taken for cluster sizes 4,5,6 where initial parameters are: {μ,σ} = {CoG, sqrt(Var(CoG)) }
 * then it does a few step least-square optimization around that. */
template<u32 N,
	typename std::enable_if<(N < FOOTClusterFit::LARGE_CLUSTER_CUTOFF)>::type* = nullptr
> 
void FitFOOTCluster(TFOOTCalProc& self, FOOTClusterFit* dest) {
	static_assert(N > 3, "Called with cluster size < 3");
	Eigen::Matrix<double,N,1> xs, ys;
	for(u32 i=0; i<N; ++i) {
		xs[i] = static_cast<double>( self._buf[i].x );
		ys[i] = self._buf[i].e;
	}
	auto params = FastGaussFit<N>(xs, ys);
	dest->a0 = params.a0;
	dest->mu = params.mu;
	dest->sigma = params.sigma;
	/* Delta assigned in parent call. */
}

/* Assumption is that the `_buf` buffer is sorted in the position already. 
 * And is also guarded that the maximum element (in energy) is index'ed as [1]. 
 * 3 points form a perfect Gauss (chi^2 undefined), and can be calculated analytically. */
template<>
void FitFOOTCluster<3>(TFOOTCalProc& self, FOOTClusterFit* dest) {
#ifdef FIT_ALL_WITH_3_POINTS
	auto it_max = std::max_element(self._buf, self._buf + self._cl_cnt, 
		[](const auto& l, const auto& r) { 
			return l.e < r.e; 
	});
	auto [x_0, e_0] = *it_max;
	int i0 = std::distance(it_max, self._buf);
	if(i0 == 0 || i0 == (int)(self._cl_cnt -1)) return;
	double e_left  = (it_max - 1)->e;
	double e_right = (it_max + 1)->e;
#else
	auto [x_0, e_0] = self._buf[1];
	double e_left  = self._buf[0].e;
	double e_right = self._buf[2].e;
#endif

	double r_left  = log( e_0 / e_left );
	double r_right = log( e_0 / e_right );
	double sigma_sq = 1 / (r_left + r_right);

	dest->sigma = sqrt(sigma_sq);

	double x = r_right / r_left;

	dest->delta = 0.5 * (1-x)/(1+x);
	dest->mu    = x_0 + dest->delta;
	dest->a0    = e_0 * exp(0.5 * dest->delta * dest->delta / sigma_sq );
}

/* For large clusters, just report the fit parameters as {μ,σ} = {CoG, sqrt(Var(CoG)) } */
template<u32 N,
	typename std::enable_if<(N == FOOTClusterFit::LARGE_CLUSTER_CUTOFF)>::type* = nullptr
> 
void FitFOOTCluster(TFOOTCalProc& self, FOOTClusterFit* dest) {
	/* Do brute force mu, sigma, a0 calculation. */
	double sy = 0.;
	double sxy = 0.;
	double sx2y = 0.;
	for(u32 i=0; i < self._cl_cnt; ++i) {
		sy  += self._buf[i].e;
		double _tmp =  self._buf[i].x * self._buf[i].e;
		sxy += _tmp;
		sx2y += self._buf[i].x * _tmp;
	}
	double mu = sxy/sy;
	double sigma_squared = std::max(sx2y/sy  - mu*mu, 0.0);
	double sigma = sqrt(sigma_squared);
	
	double phi_y = 0.;
	double phi2 = 0.;
	for(u32 i=0; i < self._cl_cnt; ++i) {
		double dx = self._buf[i].x - mu;
		double phi = exp(-0.5 * dx*dx / sigma_squared);
		phi_y += phi * self._buf[i].e;
		phi2 += phi*phi;
	}
	double a0 = phi_y / phi2;
	dest->mu    = mu;
	dest->sigma = sigma;
	dest->a0    = a0;
}

std::ostream& operator<<(std::ostream& os, const TFOOTCalProc::TClustHit& cl) {
	return os << "{" << cl.x << ", " << cl.e << "}"; 
}

void TFOOTCalProc::PrintBuff() {
	mnd_output_homogeneous_range_(std::cerr, _buf, _cl_cnt);	
	fprintf(stderr, "\n");
}
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

auto TFOOTCalProc::GetClusterType() -> ClusterType {
	using E = TFOOTCalProc::ClusterType;
	if(_ct == E::kGOOD)       return E::kGOOD;
	if(_ct == E::kMERGED)     return E::kMERGED;
	if(_ct == E::kFRAGMENTED) return E::kFRAGMENTED;
	if(_ct == E::kWAVEY)      return E::kWAVEY;
	
	return E::kUNKNOWN;
}

/* Input and output should be properly assigned a priori. */
TFOOTCalProc::TFOOTCalProc(TFOOTCalCont& out, TFOOTMapCont& in) : 
	TFOOTCalProc::Base(out, in)
{
	TFOOTMapCont& input = std::get<0>(this->in);
	if(!input.ped_s || input.ped_s->size() != N_STRIPS)
		ERROR("Input sigma array not initialized? State (null|entries): \'%s\' . Did you call TFOOTMapCont::Setup() ?\n",
				input.ped_s ? "null" : Form("%zu", input.ped_s->size()));

	if(out.FOOT_N < 0 || out.par.N < 0) 
		ERROR("Output object (TFOOTCalCont): \'%s\' uninitialized. Did you call ::Setup?", out.GetName());

	if(input.FOOT_N != out.FOOT_N)
		ERROR(EMPH2(%s) ": mapping is wrong. Deduced output FOOT%d (from setup input), but processor's input is FOOT%d ?\n", 
			_SELF_TYPE_CSTR, out.FOOT_N, input.FOOT_N);

	/* Taken from setup file. Parsed in `TFOOTCalCont::Init()` */
	X_CENTRE_THR = out.par.c_threshold;
	X_NEIGHB_THR = out.par.n_threshold;

	const auto& sigma = *input.ped_s;
	for(int i=0; i<N_STRIPS; ++i) {
		c_thr[i] = sigma[i] * X_CENTRE_THR; 
		n_thr[i] = sigma[i] * X_NEIGHB_THR; 
	}

	/* Bad/dead strips just label them with NAN'ed thresholds. */
	const std::vector<int>& bad_strips = *input.bad_strips;
	for(auto i : bad_strips) {
		c_thr[i] = BAD_STRIP_THRESHOLD; 
		n_thr[i] = BAD_STRIP_THRESHOLD; 
	}

	for(int i=0; i<N_STRIPS; ++i) {
		if(n_thr[i] < 0.1)
			ERROR("FOOT[%d -> %d], strip = %d, n_thr = %.2f too small.\n", out.FOOT_N, out.par.N, i, n_thr[i]);
		if(c_thr[i] < 0.1)
			ERROR("FOOT[%d -> %d], strip = %d, c_thr = %.2f too small.\n", out.FOOT_N, out.par.N, i, c_thr[i]);
	}

#define _ALTER_TITLE(x) \
	x->SetTitle(Form("%s : S=%.1f N=%1.f", x->GetTitle(), X_CENTRE_THR, X_NEIGHB_THR))
	_ALTER_TITLE(out.h1_cl_type);
	_ALTER_TITLE(out.h1_mult);
	_ALTER_TITLE(out.h1_X);
	_ALTER_TITLE(out.h1_dE);
	_ALTER_TITLE(out.h1_cl_sigma);
}

void TFOOTCalProc::ProcessEntry() noexcept {
	out.Clean();
	_e = &std::get<0>(in).inner().FOOTE[0]; /* Don't know if rebinding is really necessary... */

	if( std::isnan(_e[0]) ) return; /* Missing data; previous step marked it NAN. */
	
	/* Copy	the data over, since there's a write to `e`. */
	memcpy(e, _e, sizeof(e)); 
	int i = 0;
	/* Try to find a central 'seed' strip. */
	for(; i < N_STRIPS; ++i) {
		if(e[i] > c_thr[i]) {
			MakeACluster(i);
		}
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

	/* Sort by their position in the FOOT. */
	std::sort(_buf, _buf + _cl_cnt, [](const auto& l, const auto& r) { return l.x < r.x; });

	/* Check if cluster is 'merged'. Meaning that in the cluster,
	 * there exist two strips above C-threshold with 2 or more strips in-between 
	 * which are below C-threshold. */
	TClustHit* ref_hit = nullptr;
	for(u32 i=0; i < _cl_cnt; ++i) {
		const auto& [strip, adc_val] = _buf[i];
		if( adc_val > c_thr[strip] ) {
			if(ref_hit && (strip - ref_hit->x > 2)) {
				_ct_pos = ClusterType::kMERGED;
				break;
			}
			ref_hit = &_buf[i];
		}
	}

	/* Handle the wavyness in the cluster. If it is fragmented/merged, then it's wavy by default.
	 * Meaning that a cluster marked `kGOOD` has: first sequence (left-to-right) of collected strips 
	 * energy must be monotonically increasing, while second sequence must have monotonically 
	 * decreasing energy values, as to yield a proper hit structure. */
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
		out.h2_mult_e->Fill(_buf[hit].e, _cl_cnt);
	}
	
	cl_wx /= cl_e;

	/* Find the peak energy element in the cluster. */
	auto it_max = std::max_element(_buf, _buf + _cl_cnt, 
		[](const auto& l, const auto& r) { 
			return l.e < r.e; 
		});
	
	ClusterType ct = this->GetClusterType();
	out.h1_mult->Fill(_cl_cnt);
	out.h1_dE->Fill(cl_e);
	out.h1_X->Fill(cl_wx);
	out.h1_cl_type->Fill(ct);

	/* Debug cluster. */
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
	
	/* Check if it's a >3 strips cluster, then we can also do the cluster fit on the fly. */
	FOOTClusterFit cl_fit;	
	int imax = std::distance(&_buf[0], it_max); /* Index in the cluster. */
	if(imax > 0 && imax+1 < (int)_cl_cnt) {
#ifdef FIT_ALL_WITH_3_POINTS
		FitFOOTCluster<3>(*this, &cl_fit);
#else
		switch(_cl_cnt) {
			case 3: {
				FitFOOTCluster<3>(*this, &cl_fit);
				break;
			}
			case 4: {
				FitFOOTCluster<4>(*this, &cl_fit);
				break;
			}
			case 5: {
				FitFOOTCluster<5>(*this, &cl_fit);
				break;
			}
			case 6: {
				FitFOOTCluster<6>(*this, &cl_fit);
				break;
			}
			default: {
				/* Just approximate via the initial trial. Can't do much about it. */
				FitFOOTCluster<FOOTClusterFit::LARGE_CLUSTER_CUTOFF>(*this, &cl_fit);
				break;
			}
		}
		if(_cl_cnt != 3)
			cl_fit.delta = cl_fit.mu - it_max->x;
		/* Some cases cluster can generate delta outside of range [-0.5, 0.5]. 
		 * This is the ultimate null indicator field for the fit. */
		if(std::abs(cl_fit.delta) > 0.5) cl_fit.delta = NAN;
#endif
	}
	
	if(cl_fit.IsOk()) out.h1_cl_sigma->Fill(cl_fit.sigma);

	if(it_max->e > 30 and _cl_cnt == 1 and 
		out.inner()._fHeClSize1.size() == 0 and 
		out.inner()._fBadE.size() == 0) 
	{
		out.inner()._fHeClSize1.assign(e, e + N_STRIPS);
	}
	
	out.inner().AddCluster(cl_wx, cl_e, _cl_cnt, ct, cl_fit);
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
