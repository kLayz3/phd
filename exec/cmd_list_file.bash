#!/usr/bin/env bash

script_dir="$(dirname -- "$(readlink -f -- $0)")"
. "$script_dir"/common.bash
. "$script_dir"/colours.load

cmd_file=$(echo $0 | sed 's/.bash/.cmd/')

handle_h() {
	echo -e "Usage: $0 ARGS -- OPT_ARGS"
	echo -e "  -m|--make  \t\tCreate the file."
	echo -e "  -c|--cat  T\t\tCat the command list file for target T=map|cal|hit and exit."
	echo -e "  -q|--qcat T\t\tCat the command list file for target T=map|cal|hit without ❌,✅ symbols."
	echo -e "  -d|--dbg  T\t\tFor -q, -c and -e, insert /dbg/ before file name. Value should be given as bitmasks 1,2,3 "
	echo -e "             \t\tbased on if the append should go before input or output file name."
	echo -e "             \t\t1 => input, 2 => output, 3 => input+output."
	echo -e "  -e|--exec N\t\tTake Nth line of file and exec it."
	echo -e "  -n|--norun \t\tDo a dry run of the -e,--exec option."
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
do_qcat=0
do_i_dbg=0
do_o_dbg=0

do_dry=0

cat_target=""

handle_c() {
	cat_target="$1.exe"
	do_cat=1
	do_qcat=0
}
handle_q() {
	cat_target="$1.exe"
	do_cat=1
	do_qcat=1
}
handle_d() {
	case "$1" in
		[123]) ;;
		*) echo "-d : Supply a bitmask of 1, 2 or 3" >&2; exit 2;;
	esac
	do_i_dbg=$(( $1 & 1 ))
	do_o_dbg=$(( ($1 >> 1) & 1 ))
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
handle_n() {
	do_dry=1
}
handle_f() { cmd_file=$1; }

register_opt h "help"  handle_h 0
register_opt c "cat"   handle_c 1
register_opt q "qcat"  handle_q 1
register_opt d "dbg"   handle_d 1
register_opt m "make"  handle_m 0
register_opt e "exec"  handle_e 1
register_opt n "norun" handle_n 0
register_opt f "file"  handle_f 1

parse_opts "$@" || exit 2

prepare_line() {
	local line="$1" i path # 3 local variables :)
	read -r -a words <<< "$line"
	for ((i=1; i<${#words[@]}; ++i)); do
		if { [[ "${words[i-1]}" == -f ]] && (( do_i_dbg )); } ||
		   { [[ "${words[i-1]}" == -o ]] && (( do_o_dbg )); }; then
			path="${words[i]}"
			words[i]="$(dirname -- "$path")/dbg/$(basename -- "$path")"
		fi
	done
}

if (( do_cat )); then
	linenum=0
	while IFS= read -r line || [[ -n "$line" ]]; do
		((++linenum))
		prepare_line "$line"
		[[ "${words[0]}" == "./$cat_target" || "${words[0]}" == "$cat_target" ]] || continue

		for ((i=1; !do_qcat && i<${#words[@]}; ++i)); do
			if [[ "${words[i-1]}" == -f || "${words[i-1]}" == -o ]]; then
				if [[ -f "${words[i]}" ]]; then
					words[i]+=✅
				else
					words[i]+=❌
				fi
			fi
		done
		printf '%6d\t%s\n' "$linenum" "${words[*]}"
	done < "$cmd_file"
elif (( do_exec )); then
	line=$(sed -n "${linenum}p" "$cmd_file")
	prepare_line "$line"
	line="${words[*]}"
	echo -e "Replacing script with:\n$BI_YELLOW""$line" "$rest" $COLOR_OFF
	echo -e "Bye!\n"
	if(( !do_dry )); then
		exec "${words[@]}" $rest
	fi
elif (( do_create )); then
	[ -f $cmd_file ] && { echo "File: "$cmd_file" found, and will be rewritten."; ask_to_proceed; }
	
	> $cmd_file || { echo "Unable to write to "$cmd_file"; "; exit 44; }
	
	dir="${script_dir}/../merged"
	ext="main_*.root"

	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		echo "./map.exe -f merged/$fname -o mapdir/$fname -d 300" | tee -a $cmd_file
	done
	echo | tee -a $cmd_file
	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		echo "./cal.exe -f mapdir/$fname -o caldir/$fname" | tee -a $cmd_file
	done
	echo | tee -a $cmd_file	
	
	for fname in $(ls $dir/$ext); do
		fname=$(basename "$fname")
		echo "./hit.exe -f caldir/$fname -o hitdir/$fname" | tee -a $cmd_file
	done
else
	echo "Neither -c, -m, -e flag passed. Choose one."
	handle_h
	exit 22
fi
