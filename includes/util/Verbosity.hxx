#pragma once

/* ROOT libraries like Minuit, RLogger expose
 * verbosity enums but they are dependency heavy.
 * A simple variant defined here to be exposed further. 
 * Chatty variant tributed to a good friend Dr. L. Rose */

#include <optional>

enum class Verbosity : int { 
	SILENT   = 0,
	INFO     = 1,
	CHATTY   = 2,
	SPAM     = 3,
	INFINITE = 4,
};

inline bool operator==(Verbosity v, int rhs) noexcept { return static_cast<int>(v) == rhs; }
inline bool operator!=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) != rhs; }
inline bool operator< (Verbosity v, int rhs) noexcept { return static_cast<int>(v) <  rhs; }
inline bool operator> (Verbosity v, int rhs) noexcept { return static_cast<int>(v) >  rhs; }
inline bool operator<=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) <= rhs; }
inline bool operator>=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) >= rhs; }

namespace mnd {

inline std::optional<Verbosity> itov(int value) noexcept {
	switch(value) {
		case 0: return Verbosity::SILENT;
		case 1: return Verbosity::INFO;
		case 2: return Verbosity::CHATTY;
		case 3: return Verbosity::SPAM;
		case 4: return Verbosity::INFINITE;
		default: return std::nullopt;
	}
}
}
