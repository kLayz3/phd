#pragma once

#include <fstream>
#include <filesystem>
#include <type_traits>
#include <vector>
#include <string>
#include <optional>

#include <sys/wait.h>
#include <unistd.h>

namespace mnd::meta {

template<typename T>
using is_pathlike_arg =
	typename ::std::disjunction<
		std::is_base_of<std::ifstream, std::decay_t<T>>,
		std::is_constructible<std::ifstream, T&&>,
		std::is_constructible<std::string, T&&>,
		std::is_constructible<std::filesystem::path, T&&>
	>;
template<typename T>
constexpr bool is_pathlike_arg_v = is_pathlike_arg<T>::value;

} // namespace mnd::meta

namespace mnd::fs {

/* Parse a file first thru the GCC preprocessor, and then
 * try to parse the output as a sequence of lines.
 * Is not thread safe! */
auto parse_file_to_string(const std::string& ) -> std::string;
auto parse_file_to_lines (const std::string& ) -> std::vector<std::string>;

/* Make a sequence of paths based on an iterable container.
 * E.g.: path_sequence{"foo", "bar", "baz.txt"} -> "foo/bar/baz.txt" */
template<typename Range>
std::filesystem::path path_sequence(Range const& parts) {
	std::filesystem::path p;

	for(auto const& part : parts)
		p /= part;

	return p;
}

/* Resolve a path object if it's a symlink. If it is not, is just an identity function. */
std::filesystem::path resolve_maybe_symlink(const std::filesystem::path& );

template<typename T,
	typename = std::enable_if_t<
		std::is_constructible_v<std::filesystem::path, T&&>
	>
>
std::optional<std::ifstream> get_maybe_ifstream(T&& arg) {
	std::ifstream file{ std::filesystem::path{std::forward<T>(arg)} };
	if (!file.is_open())
		return std::nullopt;

	return std::optional<std::ifstream>{std::move(file)};
}

inline bool is_file_readable(const std::ifstream& file) {
	return file.is_open() && !file.fail();
}
template<typename T,
	typename = std::enable_if_t<
		std::is_constructible_v<std::filesystem::path, T&&>
	>
> bool is_file_readable(T&& arg) {
	auto file = get_maybe_ifstream(std::forward<T>(arg));
	return file && is_file_readable(*file);
}

} // namespace mnd::fs


