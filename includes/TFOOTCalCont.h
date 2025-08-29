#pragma once

#include "TContainer.h"
#include "TFOOTPedestalCont.h"

class TH1D;

class TFOOTCalCont : public TContainer {
	friend class TFOOTCalProc;
public:	
	static constexpr size_t CAPACITY = 300; //!
	static constexpr int N_STRIPS = _FOOT_N_STRIPS; //!

public:
	enum ClusterType {
		kUNKNOWN    = 0,
		kGOOD       = 1,
		kFRAGMENTED = 2, /* One or two strips missing between them. */
		kWAVEY      = 3, /* Cluster such as: _/\/\_, initial sequence isn't strictly rising, latter isn't stricly falling. */
	};
	
	enum class Orientation {
		kUNKNOWN = 0,
		kX       = 1,
		kY       = 2,
	};
	
private:
	double _x[CAPACITY] {0}; //!
	double _e[CAPACITY] {0}; //!
	double _m[CAPACITY] {0}; //!
	ClusterType _t[CAPACITY] {kUNKNOWN}; //!

public:

	Orientation o; //! 
	int FOOT_N;  //! /* Comes from sort step. */
	int POS; //! /* Force the FOOT's to be labelled 0,1,2,3,4,5,6,7 from now onward.

	TH1D* h1_mult; //!
	TH1D* h1_dE; //!
	TH1D* h2_X; //!

private:
	void AddCluster(double, double, double, ClusterType);

public:
	TFOOTCalCont();
	TFOOTCalCont(int );
	virtual ~TFOOTCalCont();

	int N;                   /* Number of clusters in the event. */
	double* fCX; //[N]       /* Cluster mean strip position. */
	double* fCE; //[N]       /* Cluster summed energy. */
	double* fCM; //[N]       /* Cluster mean multiplicity. */
	ClusterType* fCT; // [N] /* Cluster type. */

	void Clean(Option_t* option="") noexcept /* override */;
	void Init(TDictInfo info) /* override */;

	DECL_CONTAINER_METHODS
	ClassDef(TFOOTCalCont, 1);
};
