#pragma once

#include "TFOOTCalCont.h"
#include "TProcessor.h"

class TFOOTPedestalCont;

class TFOOTCalProc : public TProcessor {
public:
	static constexpr int N_STRIPS    = TFOOTCalCont::N_STRIPS; //!          /* 640 */
	static constexpr int MAX_CL_SIZE = 20;
	using ClusterType = TFOOTCalCont::ClusterType;

public:
	static constexpr double X_THR = 2.5;

	TFOOTPedestalCont* input;
	TFOOTCalCont* output;
	
	TFOOTCalProc(TFOOTPedestalCont& in, TFOOTCalCont& out) : input(&in), output(&out) {}
		
	void Init(const TDictInfo& , const TDictInfo& ) ;
	void ProcessEntry() noexcept;
			
private:
	void MakeACluster(int& );

	Double_t *e ;
	Double_t thr[N_STRIPS];
	double _buf_e[MAX_CL_SIZE] = {0};
	int _buf_i[MAX_CL_SIZE] = {0};
	int _cl_cnt ;
	
	ClusterType _ct ;

	inline void _ClearClust() noexcept { _cl_cnt = 0; _ct = ClusterType::kGOOD /* default value. */; }
	inline void _AddHit(int i) noexcept {
		_buf_e[_cl_cnt] = e[i] - thr[i];
		_buf_i[_cl_cnt] = i;
		++_cl_cnt;
	}

	bool _IsAddibleToCluster(int );
	double _CalculateMultiplicity();
};
