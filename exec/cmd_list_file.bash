#!/usr/bin/env bash
script_dir="$(dirname -- "$(readlink -f -- $0)")"
. $script_dir/common.bash

cmd_file=$(echo $0 | sed 's/.bash/.cmd/')

handle_h() {
	echo -e "Usage: $0 ARGS -- OPT_ARGS"
	echo -e "  -m|--make\t\tCreate the file."
	echo -e "  -c|--cat T\t\tCat the command list file for target T=map|cal|hit and exit."
	echo -e "  -e|--exec N\t\tTake Nth line of file and exec it."
	echo -e "\nOptional args:"
	echo -e "  -f|--file F\t\tFile name, default: "$cmd_file""
	echo -e "  OPT_ARGS\t\tPassed directly as arguments to the executed line (-e option)"
	echo
	echo "Will either take Nth line of file '$cmd_file' and exec that line"
	echo "or create the same file and exit."
	exit 0
}

do_exec=0
do_create=0
do_cat=0
cat_target=""

handle_c() {
	cat_target="$1.exe"
	do_cat=1
}
handle_m() { do_create=1; do_exec=0; }
handle_e() {
	case "$1" in
		''|*[!0-9]*) echo "-e : Supply a valid number argument "; exit 2;;
		*) linenum=$1 ;;
	esac
	do_create=0
	do_exec=1 
}
handle_f() { cmd_file=$1; }

register_opt h "help" handle_h 0
register_opt c "cat"  handle_c 1
register_opt m "make" handle_m 0
register_opt e "exec" handle_e 1
register_opt f "file" handle_f 1

parse_opts "$@"

if (( do_cat )); then
	exec cat -n $cmd_file | grep $cat_target
elif (( do_exec )); then
	line=$(sed -n "${linenum}p" "$cmd_file")
	line="${line//❌}"
	line="${line//✅}"
	echo "Replacing script with: "$line $rest""
	echo -e "Bye!\n"
	exec $line $rest
elif (( do_create )); then
	[ -f $cmd_file ] && { echo "File: "$cmd_file" found, and will be rewritten."; ask_to_proceed; } 
	
	> $cmd_file || { echo "Unable to write to "$cmd_file"; "; exit 44; }
	
	dir="${script_dir}/../merged"
	ext="main_*.root"

	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		[ -f mapdir/$fname ] && { has_output=✅; } || { has_output=❌; }
		echo "./map.exe -f merged/$fname -o mapdir/$fname$has_output -d 300" | tee -a $cmd_file
	done
	echo | tee -a $cmd_file	
	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		[ -f mapdir/$fname ] && { has_input=✅; } || { has_input=❌; }
		[ -f caldir/$fname ] && { has_output=✅; } || { has_output=❌; }
		echo "./cal.exe -f mapdir/$fname$has_input -o caldir/$fname$has_output" | tee -a $cmd_file
	done
	echo | tee -a $cmd_file	
	
	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		[ -f caldir/$fname ] && { has_input=✅; } || { has_input=❌; }
		[ -f hitdir/$fname ] && { has_output=✅; } || { has_output=❌; }
		echo "./hit.exe -f caldir/$fname$has_input -o hitdir/$fname$has_output" | tee -a $cmd_file
	done
else
	echo "Neither -c, -m, -e flag passed. Choose one."
	handle_h
	exit 22
fi
