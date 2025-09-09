#pragma once

/* Main class used for FOOT analysis and removing the pedestals. NOTE: this class is ONLY used for offline,
 * as it makes huge performance overheads to make a separate entry loop per each instance of the class.
 * This was done to help the code readability when exported to different other projects. */

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TTree.h"
#include "libs.hh"
#include "TProcessor.h"
#include "TFOOTPedestalCont.h"
#include "TOnce.hxx"
#include <utility>
#include <filesystem>
#include <fstream>
#include <type_traits>
#include "nlohmann/json.hpp"

class TFOOTPedestalProc : public TProcessor {
public:
	static constexpr int N_STRIPS          = TFOOTPedestalCont::N_STRIPS;          /* 640 */
	static constexpr int N_ASIC            = TFOOTPedestalCont::N_ASIC;            /* 10  */
	static constexpr int N_STRIPS_PER_ASIC = TFOOTPedestalCont::N_STRIPS_PER_ASIC; /* 64  */
	
	static constexpr double BAD_STRIP_CUTOFF_HI = 2.8;	
	static constexpr double BAD_STRIP_CUTOFF_LO = 1.2;
public:
	enum ProcessType {
		kGPED,
		kEPED
	} process_type = kGPED;
		
	TFOOTPedestalCont* data;

	TFOOTPedestalProc(TFOOTPedestalCont& data) :  
		data(&data) {}

	/* Rule-of-five: if a destructor or a custom move ctor/assignment
	 * is declared, then all three (dtor/move ctor/move assignment) must be declared,
	 * even if `= default;`. The copies are implicitly deleted, due to the base class. */
	 /* 
	TFOOTPedestalProc(TFOOTPedestalProc&& ) noexcept = default;
	TFOOTPedestalProc& operator=(TFOOTPedestalProc&& ) noexcept = default;
	~TFOOTPedestalProc() = default;
	*/

	void ProcessGlobalPedestal() noexcept;
	void CalcGlobalPedestal();
	void ProcessEventPedestal() noexcept;
	void CalcFinalPedestal();

	void ProcessEntry() noexcept ;

	inline Int_t Write() override { return data->Write(); }

	int GetRawADC(int chn);
	static constexpr inline int GetChannel(int asic, int offset) { return asic * N_STRIPS_PER_ASIC + offset; }

	static nlohmann::json _bad_strips; //!
	
	template <class T>
	static constexpr bool _is_accepted_arg_v =
		std::is_base_of_v<std::ifstream, std::decay_t<T>> ||
		std::is_constructible_v<std::ifstream, T&&> ||
		std::is_constructible_v<std::string, T&&> ||
		std::is_constructible_v<std::filesystem::path, T&&>;

	/**
	 * Pass either std::ifstream file handles, or filesystem paths or just raw strings
	 * Each argument represents a separate JSON file.
	 * JSON file inputs must be ASCII, not binary. Parser will throw on binary.
	 */
	template<typename T, typename... Ts>
	static void LoadBadStripsFile(T&& f0, Ts&&... other) {
		static_assert(_is_accepted_arg_v<T>,
				"LoadBadStripsFile args must be std::istream or a path-like "
				"(const char*, std::string, std::filesystem::path, std::string_view)");
		using U = std::decay_t<T>;
		namespace fs = std::filesystem;
		std::ifstream f{};

		if constexpr(std::is_base_of_v<std::istream, U>)
			f = f0;
		else if constexpr(
				std::is_constructible_v<std::ifstream, T&&> ||
				std::is_constructible_v<fs::path, U>)
			f = std::ifstream(std::forward<T>(f0));
		
		else if constexpr(std::is_constructible_v<std::string, T&&>)
			f = std::ifstream(std::string(std::forward<T>(f0)));
		
		if(!f.is_open()) 
			ERROR("Cannot open JSON file. It probably doesn't exist or unauthorized read.");

		try {
			json j = json::parse(f);
			::append_flat_json(TFOOTPedestalProc::_bad_strips, j);	
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
	std::array<double, N_STRIPS_PER_ASIC> ped_offset = {0.0};
};
