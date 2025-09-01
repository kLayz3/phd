#include "TFOOTCalProc.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>

void TFOOTCalProc::Init(const TDictInfo& info_in, const TDictInfo& info_out) {
	input->Init(info_in);
	output->Init(info_out);
	e = &input->FOOTE[0];
	if(!input->gr_s1 || input->gr_s1->GetN() != N_STRIPS)
		ERROR("Input sigma graph not initialized? State (null|entries): \'%s\'",
			input->gr_s1 ? "null" : Form("%d", input->gr_s1->GetN()));
	
	std::memcpy(thr, input->gr_s1->GetX(), sizeof(thr)); 

	std::transform(std::begin(thr), std::end(thr), std::begin(thr), [](auto x) { return TFOOTCalProc::X_THR * x; });
}

void TFOOTCalProc::ProcessEntry() noexcept {
	output->Clean();
	if( std::isnan(e[0]) ) return; /* Missing data; in previous step marked NAN by default. */
	
	int i = 0;
	for(; i < N_STRIPS; ++i) {
		if(e[i] > thr[i])
			MakeACluster(i);
	}
}

void TFOOTCalProc::MakeACluster(int& i /* starting index. */ ) {
	_ClearClust();

	/* From the loop caller, it's assumed that previous two 
	 * strips aren't part of this specific cluster. A.K.A, they're below their threshold. 
	 * And that the current strip `i` is above its respective threshold. */

	_AddHit(i);
	
	/* The moment `e[i] > thr[i]` fails, a lookahead to i+1 occurs. 
	 * If the lookahead finds a hit, e[i] is averaged out between e[i-1] and e[i+1] or stays itself, if larger. 
	 * If it fails, means both e[i] and e[i+1] fail the condition and it constitutes a 
	 * cluster end. */
	while((++i) < N_STRIPS and (e[i] > thr[i] || _IsAddibleToCluster(i)))  {
		_AddHit(i);
	}
	
	/* Handle the wavyness in the cluster. If it is fragmented, then it's wavy by default.
	 * Meaning: first sequence of collected strips energy must be monotonically increasing,
	 * while second sequence musthave monotonically decreasing energy values, as to yield a proper
	 * hit structure, if the cluster is to be marked `kGOOD`.  */
	
	if(_ct != ClusterType::kFRAGMENTED) {
		double prev_val = 0;
		enum { rising, falling } curve = rising;

		for(int hit = 0; hit < _cl_cnt; ++hit) {
			if(curve == rising and _buf_e[hit] < prev_val) {
				curve = falling;
			}
			else if(curve == falling and _buf_e[hit] > prev_val) {
				_ct = ClusterType::kWAVEY;
				break;
			}
			prev_val = _buf_e[hit];
		}
	}
	
	/* Integrate the cluster. */
	double cl_e = std::accumulate(_buf_e, _buf_e + _cl_cnt, 0.0);
	
	/* Make a weighted average for position. */
	double cl_wx = 0;
	for(int hit = 0; hit < _cl_cnt; ++hit) 
		cl_wx += _buf_i[hit] * _buf_e[hit];
	
	cl_wx /= cl_e;
	
	double cl_m = _CalculateMultiplicity();
	output->AddCluster(cl_wx, cl_e, cl_m, _ct);
}


bool TFOOTCalProc::_IsAddibleToCluster(int i /* 1, ... , 639 */) {
	/* Look ahead by one. */
	if(i < (N_STRIPS-1) and e[i+1] > thr[i]) {
		e[i] = std::max( (e[i-1] + e[i+1]) / 2, e[i] );
		_ct = ClusterType::kFRAGMENTED;
		return true; // Caller fnc adds this hit.
	}
	return false;
}

double TFOOTCalProc::_CalculateMultiplicity() {
	double cl_max = *std::max_element(_buf_e, _buf_e + _cl_cnt);
	double cl_m = 0;
	
	for(int hit = 0; hit < _cl_cnt; ++hit) 
		cl_m += _buf_e[hit] / cl_max;

	return cl_m;
}
