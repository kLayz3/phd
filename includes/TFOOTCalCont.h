#pragma once

#include "TContainer.h"
#include "TFOOTPedestalCont.h"

class TH1D;

class TFOOTCalCont : public TContainer {
	friend class TFOOTCalProc;
public:	
	static constexpr size_t INIT_CAPACITY = _FOOT_N_STRIPS_PER_ASIC; //!
	static constexpr int N_STRIPS = _FOOT_N_STRIPS; //!

public:
	enum ClusterType {
		kUNKNOWN    = 0, /* Unqualified. */
		kGOOD       = 1, /* Good cluster. Monotonically rising ADC values to peak ADC strip, then monotonically decreasing. */
		kFRAGMENTED = 2, /* One strip in noise; is missing between two sequences of the cluster. */
		kWAVEY      = 3, /* Cluster such as: _/\/\_, initial sequence isn't strictly rising, latter isn't stricly falling. */
		kMERGED     = 4, /* When two or more non-neighbouring distinct strips pass C-threshold check, and form a cluster. */
	};
	
	enum class Orientation {
		kUNKNOWN = 0,
		kX       = 1,
		kY       = 2,
	};
	
public:

	Orientation o; //! 
	int FOOT_N;  //! /* Comes from sort step. */
	int POS; //! /* Force the FOOT's to be labelled 0,1,2,3,4,5,6,7 from now onward.

	TH1I* h1_raw_mult; //!
	TH1I* h1_mult; //!
	TH1I* h1_dE; //!
	TH1I* h1_X; //!
	TH1I* h1_cl_type; //!

	TH1I* h1_dE_m1; //!
	TH1I* h1_dE_m2; //!
	TH1I* h1_dE_m3; //!
	TH1I* h1_sn_ratio; //!

private:
	void AddCluster(double, double, double, ClusterType);

public:
	TFOOTCalCont();
	virtual ~TFOOTCalCont();

	int N;                   /* Number of clusters in the event. */
	std::vector<double> fCX; /* Cluster mean strip position. */
	std::vector<double> fCE; /* Cluster summed energy. */
	std::vector<double> fCM; /* Cluster mean multiplicity. */
	std::vector<ClusterType> fCT; /* Cluster type. */

	/* Record whole event in a vector, if we find a large cluster, for some reason. */
	std::vector<double> _fBadE; /* Size will be either 0 or 640. */
	std::vector<double> _fHeClSize1; /* Size will be either 0 or 640. */

	void Clean(Option_t* option="") noexcept /* override */;
	void Init(TDictInfo info) /* override */;

	DECL_CONTAINER_METHODS
	ClassDef(TFOOTCalCont, 1);
};
