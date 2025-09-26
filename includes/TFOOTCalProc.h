#pragma once

#include "TFOOTPedestalCont.h"
#include "TFOOTCalCont.h"
#include "TProcessor.h"

class TFOOTPedestalCont;

class TFOOTCalProc : public TProcessor {
public:
	static constexpr int N_STRIPS    = TFOOTCalCont::N_STRIPS; /* 640 */
	static constexpr int MAX_CL_SIZE = 40; /* Maximal allowed cluster size. */
	static constexpr int MASSIVE_CLUSTER_CUTOFF = 20; /* After which multiplicity a cluster is called 'massive' */
	static constexpr double BAD_STRIP_FAKE_THRESHOLD = 10'000;

	static_assert(MAX_CL_SIZE > MASSIVE_CLUSTER_CUTOFF);

	using ClusterType = TFOOTCluster::ClusterType;
	using ClusterFullType = std::pair<ClusterType, ClusterType>;

	friend bool operator==(ClusterFullType&, ClusterType );
	friend bool operator!=(ClusterFullType&, ClusterType );

	enum LookAhead {
		kPOS =  1,
		kNEG = -1,
	};
	
	struct TClustHit {
		int x;    /* Position.  */
		double e; /* ADC value. */
	};

	static constexpr double X_CENTRE_THR_STATIC = 4;
	static constexpr double X_NEIGHB_THR_STATIC = 1;
	
	static_assert(X_CENTRE_THR_STATIC > X_NEIGHB_THR_STATIC, 
		"Cannot cluster correctly if seed strip threshold cutoff is smaller than neighbouring strip's threshold (AMS paper).");

	TFOOTPedestalCont* input;
	TFOOTCalCont* output;

	double X_CENTRE_THR;
	double X_NEIGHB_THR;
	
	TFOOTCalProc(TFOOTPedestalCont& in, TFOOTCalCont& out, 
			double x_seed  = X_CENTRE_THR_STATIC, 
			double x_neigh = X_NEIGHB_THR_STATIC);

	void ProcessEntry() noexcept;
	inline Int_t Write() override { return output->Write(); }

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
				input->GetName(), output->GetName(), MAX_CL_SIZE,
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
