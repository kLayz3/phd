#pragma once

#include <utility>
#include <filesystem>
#include <fstream>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TTree.h"

#include "monad/monad.hxx"
#include "nlohmann/json.hpp"

#include "TFOOTMapCont.h"
#include "TFRSGo4Cont.hxx"

struct TFOOTMapProc : TProcessor <
	TFOOTMapCont
	(TFRSGo4Cont)
> {
	static constexpr int N_STRIPS          = TFOOTMapCont::N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = TFOOTMapCont::N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = TFOOTMapCont::N_STRIPS_PER_ASIC; /* 64  */
	
	static constexpr double BAD_STRIP_CUTOFF_HI = 2.8;	
	static constexpr double BAD_STRIP_CUTOFF_LO = 1.2;

	static constexpr int N_BATCH_PEDESTAL = 20000;

	using Base = TProcessor<TFOOTMapCont(TFRSGo4Cont)>;

public:
	enum ProcessType {
		kINITIAL_BATCH,
		kMAIN
	} process_type = kINITIAL_BATCH;
	enum class CableSwapped { 
		NO = 0, 
		YES = 1
	} is_swapped = CableSwapped::NO;
	
	struct NBatchPedestal {
		u32 amount;
	} n_batch_pedestal;
	
	double dt_veto = NAN; /* Optionally given. */

	int N; // Corresponds to FOOT_N from the output container.
	uint32_t nsampled = 0;

	TFOOTMapProc() = default;
	TFOOTMapProc(TFOOTMapCont& out, const TFRSGo4Cont& in, 
		NBatchPedestal n_batch = { N_BATCH_PEDESTAL },
		CableSwapped is_cabling_swapped = CableSwapped::NO, double dt_veto_ = NAN);

	void ProcessInitialPedestal() noexcept;
	
	void ProcessEventPedestal() noexcept;
	void CalcGlobalPedestal();
	void CalcFinalPedestal();

	void ProcessEntry() noexcept;

	int GetRawADC(int chn);
	static constexpr inline int GetChannel(int asic, int offset) { return asic * N_STRIPS_PER_ASIC + offset; }

	static nlohmann::json _bad_strips; //!
	
	/**
	 * Pass either std::ifstream file handles, or filesystem paths or just raw strings
	 * Each argument represents a separate JSON file.
	 * JSON file inputs must be ASCII, not binary. Parser will throw on binary.
	 */
	template<typename T, typename... Ts>
	static void LoadBadStripsFile(T&& arg, Ts&&... other) {
		static_assert(mnd::is_pathlike_arg_v<T>,
				"LoadBadStripsFile args must be std::ifstream or a path-like "
				"(const char*, std::string, std::filesystem::path, std::string_view)");

		std::ifstream f{}; 
		try {
			f = mnd::get_maybe_ifstream(std::forward<T>(arg)).value();
		} catch(const std::exception& e) {
			ERROR("Cannot open JSON file. It probably doesn't exist or unauthorized read.");
		}

		try {
			nlohmann::json j = nlohmann::json::parse(f);
			append_flat_json(TFOOTMapProc::_bad_strips, j);	
		} catch(std::exception const& e) {
			WARN("Json parsing failed. Reason: %s\n", e.what());
			throw;
		}
		(LoadBadStripsFile(other), ...);
	} 

	/**
	 * Call this at the end of the processing, to also write to the `bad_strips`
	 * the strips passed by the user, which are stored in the static JSON field.
	 */
	int ParseStaticBadStrips();

private:
	std::array<double, N_STRIPS> initial_gped = {0.0};
	std::array<double, N_STRIPS> current_gped = {0.0};
	
	bool initial_calculated = false;
	u32 batch_index = 0;

	std::array<double, N_STRIPS_PER_ASIC> ped_offset = {0.0};

	struct FOOTView {
		bool TSBad; u32 TLO; u32 THI; u32 SY;
		u32 Ndata; u32* E;
	};
	Scaler<32> init_timing {};

	static const FOOTView GetPtrs(const TFRSSortEvent* e, int N);

	static void parse_json_string(std::vector<int>& out, std::string s);
	static void parse_json_as_int_vec(std::vector<int>& out, const nlohmann::json& j);
	static void append_flat_json(nlohmann::json& dst, const nlohmann::json& src);
};
