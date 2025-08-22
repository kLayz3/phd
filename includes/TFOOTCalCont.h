#pragma once

#include "TContainer.h"
#include "TFOOTPedestalCont.h"

class TFOOTCalCont : public TContainer {
public:	
	static constexpr int N_STRIPS = _FOOT_N_STRIPS; //!
	static constexpr int N_ASIC = _FOOT_N_ASIC; //!
	static constexpr int N_STRIPS_PER_ASIC = _FOOT_N_STRIPS_PER_ASIC; //!

	static constexpr size_t CAPACITY = 100; //!
	
public:
	enum ClusterType {
		kUNKNOWN    = 0,
		kGOOD       = 1,
		kFRAGMENTED = 2, /* One or two strips missing between them. */
		kWAVEY      = 3, /* Cluster such as: _/\/\_ */
	};
	
	enum class Orientation {
		kUNKNOWN = 0,
		kX       = 1,
		kY       = 2,
	};
private:
	double _x[CAPACITY] {0}; //!
	double _y[CAPACITY] {0}; //!
	double _m[CAPACITY] {0}; //!
	ClusterType _c[CAPACITY] {kUNKNOWN}; //!
public:
	Orientation o; //!
 
	int FOOT_N;  //! /* Comes from sort step. */
	int FOOT_ID; //! /* Force the FOOT's to be called X
	TFOOTCalCont();
	TFOOTCalCont(int );
	virtual ~TFOOTCalCont();

	int N; 
	double* fCx; //[N]
	double* fCy; //[N]
	u32* fCM; //[N]
	ClusterType* fCt; // [N]

	ClassDef(TFOOTCalCont, 1);
};
