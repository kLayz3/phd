#!/usr/bin/env -S fish --no-config

set -l usage "Usage: $(basename (status filename)) INPUT.root [OUTPUT.pdf]"
set -l help "
Delegates to the `draw_transport.C` macro, which will just draw
(and optionally generate a PDF file) from the current 'bare' transport process.
E.g. just shoots particles and measures the scattering statuses after each layer of FOOT."

if contains -- -h $argv; or contains -- --help $argv
	printf '%s\n%s\n' "$usage" "$help"
	exit 0
end

set -l narg (count $argv)
if test $narg -gt 2; or test $narg -lt 1
	printf '%s\n' "$usage" >&2
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
	root -l "$script_dir/transport/draw_transport.C($call_args)"
