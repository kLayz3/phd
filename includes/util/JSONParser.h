#pragma once

#include "json_struct_def.hh"

/** 
 * Parse a file first thru the GCC preprocessor, and then
 * try to parse the output as json. Returns a variant of the
 * valid JSON or a corresponding error code. Is not thread safe!
 */
nlohmann::json ParseJSON(const std::string& );
