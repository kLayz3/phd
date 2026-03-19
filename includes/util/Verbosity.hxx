/* ROOT libraries like Minuit, RLogger expose
 * verbosity enums but they are dependency heavy.
 * A simple variant defined here to be exposed further. 
 * Chatty variant tributed to a good friend Dr. L. Rose */

#pragma once

enum class Verbosity : int { 
	SILENT = 0,
	INFO   = 1,
	CHATTY = 2
};

inline bool operator==(Verbosity v, int rhs) noexcept { return static_cast<int>(v) == rhs; }
inline bool operator!=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) != rhs; }
inline bool operator< (Verbosity v, int rhs) noexcept { return static_cast<int>(v) <  rhs; }
inline bool operator> (Verbosity v, int rhs) noexcept { return static_cast<int>(v) >  rhs; }
inline bool operator<=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) <= rhs; }
inline bool operator>=(Verbosity v, int rhs) noexcept { return static_cast<int>(v) >= rhs; }
