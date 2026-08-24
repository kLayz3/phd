#pragma once

#include "monad/monad.hxx"
#include "TFRSMapCont.h"
#include "TFRSCalCont.h"

//#define TFRSCALPROC_VERBOSE_
#define TFRSCALPROC_SINGLEHIT

struct TFRSCalProc : TProcessor <
	TFRSCalCont
	(TFRSMapCont, TTrigMapCont)
> {
	using Base = TProcessor<TFRSCalCont(TFRSMapCont, TTrigMapCont)>;
	
	TFRSCalProc(TFRSCalCont& , const TFRSMapCont&, const TTrigMapCont& );
	TFRSCalProc() = default;

#ifdef TFRSCALPROC_VERBOSE_
	inline static u64 called_ = 0;
	inline static u64 errors_ = 0;
	~TFRSCalProc() { WARN("N_ERRORS: " EMPH1(%lu\n), errors_); }
#endif

	constexpr static auto N_VALID_TPC = RNFRSCal::N_VALID_TPC;
	constexpr static auto N_VALID_SCI = RNFRSCal::N_VALID_SCI;

	void ProcessEntry() noexcept;

private:
	/* Encapsulating viable data from single TPC, single anode channel */
	struct TPCHitCandidate {
		i32 a_tdc;
		i32 dl_tdc;
		i32 dr_tdc;
		enum class Status { kFINE, kIN_CONFLICT, kDROPPED } status = Status::kFINE;

		TPCHitCandidate() = default;
		TPCHitCandidate(i32 _a_tdc, i32 _dl_tdc, i32 _dr_tdc, i32 _ref_tdc) noexcept :
			a_tdc(_a_tdc), dl_tdc(_dl_tdc), dr_tdc(_dr_tdc) {}
		
		inline i32 CSum() { return dl_tdc + dr_tdc - 2*a_tdc; }
	};
	using TPCHitCandidateList = std::vector<TPCHitCandidate>;
	constexpr static std::size_t CANDIDATE_LIST_CAPACITY = RNTPCMap::MAX_SIZE * RNTPCMap::MAX_SIZE * RNTPCMap::MAX_SIZE;
	
	using TPCConflict = std::pair<TPCHitCandidate*, TPCHitCandidate*>;
	using TPCConflicts = std::vector<std::pair<TPCHitCandidate*, TPCHitCandidate*>>;

	bool TryResolveViaAnodeDiff(TPCConflict&, std::array<double, 2>& limits);

	enum class MeasurementType { kLEFT, kRIGHT, kANODE };
	template<MeasurementType m> 
	bool TryResolveViaDiff(TPCConflict&, const std::array<double, 2>& limits, const TPCHitCandidateList& other, bool );

	struct TPCHitCandidateExtended {
		std::array<i32,2> a_tdc;
		i32 dl_tdc;
		i32 dr_tdc;
		i32 ref_tdc;
		TPCHitCandidate hit; 
		TPCHitCandidateExtended() = default;
		TPCHitCandidateExtended(i32 _a_tdc0, i32 _a_tdc1, i32 _dl_tdc, i32 _dr_tdc, i32 _ref_tdc) noexcept :
			a_tdc{_a_tdc0, _a_tdc1}, dl_tdc(_dl_tdc), dr_tdc(_dr_tdc), ref_tdc(_ref_tdc) {}

		inline bool operator<(const TPCHitCandidateExtended& rhs) const noexcept {
			return ref_tdc < rhs.ref_tdc;
		}
	};
	using TPCHitCandidateExtendedList = std::vector<TPCHitCandidateExtended>;
	
	/* One list per anode, per delay-line. */
	std::array<
		std::array<TPCHitCandidateList, 2>, 2
	> candidate_list {};

	/* One list per anode, per delay-line. */
	std::array<
		std::array<TPCConflicts, 2>, 2
	> conflicts {};

	/* One list per delay-line in a TPC. */
	std::array<TPCHitCandidateExtendedList, 2> full_candidate_list {};

	void ProcessTPC(int ) noexcept;
	void PreProcessTPC(int ) noexcept;
	void ProcessCSum(int, int ) noexcept;
	void ResolveTPCDelayLineConflicts(int, int ) noexcept;
	void FilterRemainingConflicts() noexcept;
	void ProcessTPC_YCut(int ) noexcept;
	void PostProcessTPC(int ) noexcept;

#ifdef TFRSCALPROC_SINGLEHIT
	void ProcessSingleHit(int, int ) noexcept;
#endif

	/* === SCI analysis helper fnc's. === */
	void ProcessSci(int ) noexcept;
};
