#pragma once

#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include "monad/monad.hxx"

struct TFOOTMapCont;

struct TFOOTCalProc : TProcessor <
	TFOOTCalCont
	(TFOOTMapCont)
> {
	using Base = TProcessor<TFOOTCalCont(TFOOTMapCont)>;
	static constexpr int N_STRIPS    = TFOOTCalCont::N_STRIPS; /* 640 */
	static constexpr int MAX_CL_SIZE = 40; /* Maximal allowed cluster size. */
	static constexpr int MASSIVE_CLUSTER_CUTOFF = 20; /* After which multiplicity a cluster is called 'massive' */
	static constexpr double BAD_STRIP_FAKE_THRESHOLD = NAN;

	static_assert(MAX_CL_SIZE > MASSIVE_CLUSTER_CUTOFF);

	using ClusterType = RNFOOTCluster::ClusterType;
	using ClusterFullType = std::pair<ClusterType, ClusterType>;

	friend bool operator==(ClusterFullType&, ClusterType );
	friend bool operator!=(ClusterFullType&, ClusterType );

	enum LookAhead { kNEG, kPOS, };
	
	/* Single hit belonging to a cluster. */
	struct TClustHit {
		int x;    /* Position.  */
		double e; /* ADC value. */
	};

	double X_CENTRE_THR;
	double X_NEIGHB_THR;

	TFOOTCalProc(TFOOTCalCont& out, TFOOTMapCont& in);
	TFOOTCalProc() = default;

	void ProcessEntry() noexcept;

	template<LookAhead> bool _IsAddibleToCluster(const int );

private:
	void MakeACluster(int& );

	Double_t* _e;
	Double_t e[N_STRIPS]{0};
	Double_t c_thr[N_STRIPS] = {0};
	Double_t n_thr[N_STRIPS] = {0};

	/* Temporarily storing the information
	 * for a specific collected cluster. */

	TClustHit _buf[MAX_CL_SIZE]{};

	u32 _cl_cnt = 0; /* Points to last valid index in `_buf`.  */	
	ClusterFullType _ct { ClusterType::kUNKNOWN, ClusterType::kUNKNOWN };
	ClusterType GetClusterType();

	inline void _ClearClust() noexcept { _cl_cnt = 0; _ct = {ClusterType::kGOOD, ClusterType::kGOOD}; }
	inline void _AddHit(int i) noexcept {
		if(_cl_cnt >= MAX_CL_SIZE) {
			WARN("(%s => %s) cluster size exceeded: " EMPH(max = %d) ". Ignoring the hit (strip:ADC): " EMPH(%d: %7.2f; N-Thr: %5.2f (X=%.1f)) ". Consider increasing the buffer size.\n", 
				std::get<0>(in).GetName(), out.GetName(), MAX_CL_SIZE,
				i, e[i], n_thr[i], X_NEIGHB_THR);
			return;
		}
		_buf[_cl_cnt].x = i;
		_buf[_cl_cnt].e = e[i];
		++_cl_cnt;
	}

	/* Save the state of where the possible cluster fragmentation took place. */
	TClustHit _frag{};
};
