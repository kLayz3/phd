#pragma once

#include <charconv>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <type_traits>
#include <vector>

#include "monad/accessories/typedefs.hxx"
#include "util/Option.hxx"

/* Header that organises all the text/string parsing. */

/* File names are often of the form: `main_0XXX_0YYY.root`, as such
 * metadata'ing multiple files can be concatenated e.g.:
 * => main_0123_0144.root
 * +  main_0145_0174.root
 * +  main_0175_0178.root
 * ----------------------
 * =  main_0123_0178.root
 * */

namespace mnd {

using namespace std::literals;

inline constexpr const char* file_prefix = "main";
inline constexpr u32 nchars_run_number = 4;

/* Generally, we expect files to be called with this PCRE2 format. E.g.:
 * main_0051_0085_dbg.root
 * Group 1: 0–4   main
 * Group 2: 5–9   0051
 * Group 3: 10–14 0085 */
inline constexpr const char* rootfile_pattern
	= R"(^(\w.+?)_(\d{1,4})(?:_(\d{1,4}))?(?:.+?)?\.root$)";

/* If "main_0023_0144.root" => "main"s */
std::string run_name(std::string_view );

/* If "main_0023_0144.root" => 23 */
u32 rootfile_start_number(std::string_view );

/* If "main_0023_0144.root" => 144 */
u32 rootfile_end_number(std::string_view );

/* If "main_0023_0144.root" => { 23, 144 } */
std::pair<u32, u32> rootfile_number_bounds(std::string_view );

/* If "main_0023_0144.root" => { "main"s, 23, 144 } */
std::tuple<std::string, u32, u32>
rootfile_info(std::string_view );

/* If { "main_0023_0144.root", "main_0145_0174.root" } => { 23, 174 }.
 * It doesn't internally sort the sequence. Assumes sequence comes already sorted. */
std::pair<u32, u32> rootfile_number_bounds(const std::vector<std::string>& );

/* If `{ main_0023_0144.root, main_0145_0174.root }` => "main_0023_0174". Doesn't have the extension.
 * It doesn't internally sort the sequence. Assumes sequence comes already sorted. */
std::string rootfile_names_concatenated(const std::vector<std::string>& );

/* One liner to return a stem of the string representing a possible path.
 * If "caldir/foo/main_0012.root"sv => "main_0012" */
std::string file_stem(std::string_view );

/* {"hehe.txt", "heha.jpg", "hehe.mk"} -> "heh"s */
auto longest_common_prefix(::mnd::span<const std::string> ) -> std::string;

/* Check if the text (1st arg) has a specific prefix (2nd arg).
 * If the prefix arg is sized 0, returns false (stupid check). */
auto has_prefix(std::string_view , std::string_view ) -> bool;

/* Split string according to the separator, whilest minding
 * single/double quotes literal text. Also trimming the whitespace on both sides,
 * and treating backslash as escaping the next character everywhere.
 * A literal backslash must be doubled; a dangling backslash throws.
 * Quoted or escaped whitespace is preserved.
 * Example: split( R"(Bleebloop\,,   '1,2,3'"4,5", '(1,2,3)\'')"sv )
 * => ["Bleebloop,", "1,2,34,5", "(1,2,3)'"] */
[[nodiscard]] auto split(std::string_view , char ) -> std::vector<std::string>  ;

/* Same rules as split(), just w/ views instead. Passing rvalue here will dangle badly.
 * Since this cannot mutate the original view, the escape chars are still present here.
 * Example: split_view( R"(Bleebloop\,,   '1,2,3','(1,2,3)\'')" )
 * => [R"(Bleebloop\,)"sv, R"(   '1,2,3')"sv, R"('(1,2,3)\'')"sv] */
[[nodiscard]] auto split_view(std::string_view , char ) -> std::vector<std::string_view>;

/* Return a vector of views to the underlying sequence of strings. */
auto to_views(const std::vector<std::string>& ) -> std::vector<std::string_view> ;

/* Remove leading and trailing whitespaces from a view. Doesn't mutate the underlying buffer,
 * only returns the trimmed view object. */
inline constexpr std::string_view trim(std::string_view sv) noexcept {
	constexpr auto whitespace = " \t\n\r\f\v";

	const auto first = sv.find_first_not_of(whitespace);
	if(first == std::string_view::npos)
		return {};

	const auto last = sv.find_last_not_of(whitespace);
	return sv.substr(first, last - first + 1);
}

} // namespace mnd

/* Small DSL to properly tokenize and parse SECTION(name) { ... } blocks
 * used to pass command-line args inside programs. */

namespace mnd::sdsl {

/* Parse top-level SECTION(label) { body } blocks from the text. */
auto extract_text_body(
	std::string_view , // SECTION
	std::string_view , // section name
	std::string_view   // full text
) -> mnd::Option<std::string_view>;

/* From the text, exact all SECTION(label) . */
auto extract_available_sections(
	std::string_view , // SECTION
	std::string_view   // full text
) -> std::vector<std::string_view>;

} // namespace mnd::sdsl
