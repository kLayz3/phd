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
		i32 a_tdc;
		i32 dl_tdc;
		i32 dr_tdc;
		i32 ref_tdc;
		
		TPCHitCandidate() = default;
		TPCHitCandidate(i32 _a_tdc, i32 _dl_tdc, i32 _dr_tdc, i32 _ref_tdc) noexcept :
			a_tdc(_a_tdc), dl_tdc(_dl_tdc), dr_tdc(_dr_tdc), ref_tdc(_ref_tdc) {}
	};

	struct TPCHitCandidateExtended {
		std::array<i32,2> a_tdc;
		i32 dl_tdc;
		i32 dr_tdc;
		i32 ref_tdc;
		TPCHitCandidate hit; 
		TPCHitCandidateExtended() = default;
		TPCHitCandidateExtended(std::array<i32,2> _a_tdc, i32 _dl_tdc, i32 _dr_tdc, i32 _ref_tdc) noexcept :
			a_tdc(std::move(_a_tdc)), dl_tdc(_dl_tdc), dr_tdc(_dr_tdc), ref_tdc(_ref_tdc) {}

		inline bool operator<(const TPCHitCandidateExtended& rhs) const noexcept {
			return ref_tdc < rhs.ref_tdc;
		}
	};

	constexpr static std::size_t CANDIDATE_LIST_CAPACITY = RNTPCMap::MAX_SIZE * RNTPCMap::MAX_SIZE * RNTPCMap::MAX_SIZE;
	using TPCHitCandidateList         = std::vector<TPCHitCandidate>;
	using TPCHitCandidateExtendedList = std::vector<TPCHitCandidateExtended>;
	
	/* One list per anode in a delay-line. */
	std::array<TPCHitCandidateList, 2>         candidate_list {};
	/* One list per delay-line in a TPC. */
	std::array<TPCHitCandidateExtendedList, 2> full_candidate_list {};

	void PreProcessTPC(int ) noexcept;
	void ProcessDelayLine(int, int ) noexcept;
	void PostProcessTPC(int ) noexcept;

	void ProcessS2Angle() noexcept;
	void ProcessS4Angle() noexcept;

	/* === SCI analysis helper fnc's. === */
	void ProcessSci(int ) noexcept;
};
