#!/usr/bin/env -S fish --no-config

set -l usage "Usage: $(basename (status filename)) INPUT.root"
set -l help "
Delegates to the `print_final_status.C` macro, which will report
the status of the currently 'guessed' resolutions in the analysis simulation,
and compare it to the real measured resolutions"

if contains -- -h $argv; or contains -- --help $argv
	printf '%s\n%s\n' "$usage" "$help"
	exit 0
end

if test (count $argv) != 1
	printf '%s\n%s\n' "$usage" "$help" >&2
	exit 2
end

set -l script_dir (dirname (realpath (status -f)))
set -l cpp_args
for arg in $argv
	# Escape each argument for a C++ string literal, without eval.
	set -l escaped (string replace -a -- '\\' '\\\\' "$arg")
	set escaped (string replace -a -- '"' '\\"' "$escaped")
	set -a cpp_args "\"$escaped\""
end

set -l call_args ""
if test (count $cpp_args) -gt 0
	set call_args (string join ', ' -- $cpp_args)
end

exec env ROOT_INCLUDE_PATH="$script_dir/../includes:$ROOT_INCLUDE_PATH" \
	root -l -b -q "$script_dir/transport/print_final_status.C($call_args)"
