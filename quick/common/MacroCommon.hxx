#pragma once

#include "util/MacroHelpers.h"
#include "TFRSCalCont.h"

struct TPCRef {
    static constexpr u32 N_UPSTREAM_TPC   = TPCParam::N_UPSTREAM_TPC; // 3
    static constexpr u32 N_DOWNSTREAM_TPC = TPCParam::N_DOWNSTREAM_TPC; // 1
    static constexpr u32 N_S2_TPC = N_UPSTREAM_TPC + N_DOWNSTREAM_TPC;
    static constexpr auto TPC_WIDTH = TPCParam::TPC_WIDTH;

    /* Returns true on a good setup. */
	bool SetN(u32 x) {
		if(x < N_S2_TPC) { n = x; return true; }
		const std::string xs = std::to_string(x);
		try {
			u32 i = RNFRSCal::tpc_moniker.at( xs.c_str() );
			if(i <= N_UPSTREAM_TPC) { 
                n = i; 
                return true; 
            }
		}
		catch(...) {
            WARN("Bad query for TPC moniker, key: \'%s\'\n", xs.c_str());
        }
		return false;
	}
	operator bool() const { return n < N_S2_TPC && (use[0] || use[1]); } 
    bool IsUpstream() const { return n < N_UPSTREAM_TPC; }
    bool IsDownstream() const { return n >= N_UPSTREAM_TPC && n < N_S2_TPC; }

	u32 n; // 0 => TPC21; 1 => TPC22; 2 => TPC23;  3 => TPC24
	std::array<bool,2> use; // Which delay lines to use
};

inline std::istream& operator>>(std::istream& in, TPCRef& out) {
	/* Read an integer first */
	u32 N; unsigned char c;
	if(!(in >> N) || !out.SetN(N)) {
        WARN("Bad parse for TPCRef index number.\n");
		in.setstate(std::ios::failbit);
		return in;
	}
	
	/* Read the separator */
	if(!(in >> c) || c != ':') {
        WARN("Bad parse for TPCRef separator ':'.\n");
		in.setstate(std::ios::failbit);
		return in;
	}

	/* Read two booleans next with separator ',' */
	return ::mnd::template operator>> <','>(in, out.use);
}
inline std::ostream& operator<<(std::ostream& os, const TPCRef& out) {
	return os << '[' << out.n << out.use << ']';
}
