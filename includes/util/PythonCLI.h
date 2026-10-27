#pragma once

#include "PrettyHisto.h"

#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/* API helper string, defined in PythonCLI.cxx. */
extern const char* helper_;

namespace mnd::plot {

/* What is enough to know, for a single Pyth image, is:
 * [1] The sequence of names,
 *   which CAN be defined afterwards (just cached here, during arg parsing)
 * [2] Sequence of individual styles (same or smaller size) 
 * [3] General figure style. */

struct PyInfo {
	struct ObjectInfo {
		std::string name;
		std::vector<std::pair<std::string, std::string>> options;

		void add_object_cfg(std::string_view );
	};
	
	std::vector<ObjectInfo> objects {};
	FigureStyle fig_style {};

	void add_figure_cfg(std::string_view );
	static PyInfo parse(std::string_view );
};

// Needed for CLI11.
std::istream& operator>>(std::istream& , PyInfo& );
std::ostream& operator<<(std::ostream& , const PyInfo& );

std::ostream& operator<<(std::ostream& , const PyInfo::ObjectInfo& );

} // namespace mnd::plot

