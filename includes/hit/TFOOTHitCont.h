#include "monad/monad.hxx"
#include "TFOOTCalCont.h"

struct FOOTHit {
	double x,y,z,e;
	
	virtual ~FOOTHit();
	ClassDef(FOOTHit, 1);
};


struct RNFOOTHit {
	virtual ~RNFOOTHit();
	ClassDef(RNFOOTHit, 1);
};

struct TFOOTHitCont : TContainer<RNFOOTHit> {
	constexpr static int N_STRIPS = RNFOOTCal::N_STRIPS;
	
	TFOOTHitCont();
	void Setup() override;
	void Init(TDictInfo ) override;
};
