#pragma once

#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "core/TProcessor.hxx"

struct TFRSCalProc : TProcessor <
	TFRSCalCont
	(TFRSMapCont)
> {
	using Base = TProcessor<TFRSCalCont(TFRSMapCont)>;

	constexpr static auto N_VALID_TPC = RNFRSCal::N_VALID_TPC;
	static inline u64 ncalled = 0;
	
	TFRSCalProc(TFRSCalCont& , const TFRSMapCont& );
	TFRSCalProc() = default;

	void ProcessEntry() noexcept;

private:
	/* Encapsulating viable data from single TPC, single anode channel */
	struct TPCHitCandidate {
		int a_tdc;
		int dl_tdc;
		int dr_tdc;
		int ref_tdc;
		TPCHitCandidate(int a_tdc, int dl_tdc, int dr_tdc, int ref_tdc) noexcept :
			a_tdc(a_tdc), dl_tdc(dl_tdc), dr_tdc(dr_tdc), ref_tdc(ref_tdc) {}
		inline int CSum() noexcept { return dl_tdc + dr_tdc - (a_tdc << 1); }
	};
	struct TPCHit {
		int anode_i;
		double x;
		double y;
		TPCHit(int i, double x, double y) noexcept :
			anode_i(i), x(x), y(y) {}
	};

	constexpr static std::size_t CANDIDATE_LIST_CAPACITY = 8;
	using TPCHitCandidateList = std::vector<TPCHitCandidate>;
	
	std::array<TPCHitCandidateList, 4> candidate_list {};
	std::array<TPCHitCandidateList, 4> initial_candidate_list {};

	constexpr static std::size_t HIT_LIST_CAPACITY = 10; 
	using TPCHitList = std::vector<TPCHit>;
	TPCHitList tpc_hits;

	/* Keep a copy locally on the stack of *this* instance. 
	 * Indeed, the param array should exist here in the Proc, but we keep 
	 * original anyway in Cont to hold it in a ROOT file. */
	std::array<TPCParam, N_VALID_TPC> tpc_param;

	void ProcessTPC(int ) noexcept;
	bool IsUniqueTPCMeasurement(const TPCHitCandidateList&, const TPCHitCandidate& ) noexcept;
};
