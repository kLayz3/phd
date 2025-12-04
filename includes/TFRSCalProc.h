#pragma once

#include "TFRSMapCont.h"
#include "TFRSCalCont.h"
#include "monad/monad.hxx"

struct TFRSCalProc : TProcessor <
	TFRSCalCont
	(TFRSMapCont)
> {
	using Base = TProcessor<TFRSCalCont(TFRSMapCont)>;
	
	TFRSCalProc(TFRSCalCont& , const TFRSMapCont& );
	TFRSCalProc() = default;

	constexpr static auto N_VALID_TPC = RNFRSCal::N_VALID_TPC;
	constexpr static auto N_VALID_SCI = RNFRSCal::N_VALID_SCI;

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
		inline bool operator<(const TPCHitCandidate& rhs) const noexcept {
			return ref_tdc < rhs.ref_tdc;
		}
	};

	constexpr static std::size_t CANDIDATE_LIST_CAPACITY = 8;
	using TPCHitCandidateList = std::vector<TPCHitCandidate>;
	
	std::array<TPCHitCandidateList, 4> candidate_list {};
	std::array<TPCHitCandidateList, 4> initial_candidate_list {};

	constexpr static std::size_t HIT_LIST_CAPACITY = 10; 

	void ProcessTPC(int ) noexcept;
	bool IsUniqueTPCMeasurement(const TPCHitCandidateList&, const TPCHitCandidate& ) noexcept;

	/* === SCI analysis helper fnc's. === */
	void ProcessSci(int ) noexcept;
};
