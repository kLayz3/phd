#!/usr/bin/env bash
script_dir="$(dirname -- "$(readlink -f -- $0)")"
. $script_dir/common.bash

cmd_file=$(echo $0 | sed 's/.bash/.cmd/')

handle_h() {
	echo -e "Usage: $0 ARGS"
	echo -e "  -m|--make\t\tCreate the file."
	echo -e "  -c|--cat\t\tCat the command list file and exit."
	echo -e "  -e|--exec N\t\tTake Nth line of file and exec it."
	echo -e "\nOptional args:"
	echo -e "  -f|--file F\t\tFile name, default: "$cmd_file""
	echo
	echo "Will either take Nth line of file '$cmd_file' and exec that line"
	echo "or create the same file and exit."
	exit 0
}

do_exec=0
do_create=0
do_cat=0

handle_c() { do_cat=1; }
handle_m() { do_create=1; do_exec=0; }
handle_e() { linenum=$1; do_create=0; do_exec=1; }
handle_f() { cmd_file=$1; }

register_opt h "help" handle_h 0
register_opt c "cat"  handle_c 0
register_opt m "make" handle_m 0
register_opt e "exec" handle_e 1
register_opt f "file" handle_f 1

parse_opts "$@"
if (( !do_exec && !do_create && !do_cat )); then
	echo "Neither -c, -m, -e flag passed. Choose one."
	handle_h
	exit 22;
fi

if (( do_cat )); then
	exec cat -n $cmd_file
elif (( do_exec )); then
	line=$(sed -n "${linenum}p" "$cmd_file")
	echo "Replacing script with: "$line""
	echo -e "Bye!\n"
	exec $line
else
	[ -f $cmd_file ] && { echo "File: "$cmd_file" found, and will be rewritten."; ask_to_proceed; } 
	
	> $cmd_file || { echo "Unable to write to "$cmd_file"; "; exit 44; }
	
	dir="${script_dir}/../merged"
	ext="main_*.root"

	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		echo "./map.exe -file merged/$fname -output mapdir/$fname -foot-dt 300" | tee -a $cmd_file
	done
	echo | tee -a $cmd_file	
	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		echo "./cal.exe -file mapdir/$fname -output caldir/$fname -gm" | tee -a $cmd_file
	done
fi
