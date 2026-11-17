#!/usr/bin/env -S fish --no-config

if test (count $argv) -gt 2
	echo 'Usage: draw_transport.fish [INPUT.root [OUTPUT.pdf]]' >&2
	exit 2
end

set -l script_dir (dirname (realpath (status filename)))
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
