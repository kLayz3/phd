# Few hacks to basically split cmdline parsing in bash to multiple
# various callbacks that get registered and then later invoked.
# Check examples in the directory how to use. --k41y23

declare -A OPT_HANDLERS=()
declare -A OPT_HAS_ARG=()

SHORT_SPEC=""
LONG_SPEC=()

register_opt() {
	[ $# -ne 4 ] && { echo "[[INTERNAL]] register_opt: "$@" ; passed $# != 4 arguments."; exit 67; }

	local short="$1"
	local long="$2"
	local handler="$3"
	local has_arg="$4" # either 0 or 1, if the option requires value or not.

	OPT_HANDLERS["$short"]="$handler"
	OPT_HANDLERS["$long"]="$handler"
	OPT_HAS_ARG["$short"]="$has_arg"
	OPT_HAS_ARG["$long"]="$has_arg"

	if (( has_arg )); then
		SHORT_SPEC+="${short}:"
		LONG_SPEC+=("${long}:")
	else
		SHORT_SPEC+="${short}"
		LONG_SPEC+=("${long}")
	fi
}

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
parse_opts() {
	local parsed	
	parsed=$(getopt -o "$SHORT_SPEC" --long "$(IFS=,; echo "${LONG_SPEC[*]}")" -- "$@") || return 1
	eval set -- "$parsed"

	while true; do
		case "$1" in
			--)
				shift
				break
				;;
			-*)
				local opt="${1#-}"
				opt="${opt#-}"
				# Strip twice due to long option
				local handler="${OPT_HANDLERS[$opt]}"
				local has_arg="${OPT_HAS_ARG[$opt]}"
				if (( has_arg )); then
					"$handler" "$2"
					shift 2
				else
					"$handler"
					shift
				fi
				;;
			*)
				break
				;;
		esac
	done
}

# This function will just return a 0 on `no`. Won't promptly exit.
ask_to_proceed () {
	do_it=
	read -p "Do you wish to proceed? [yY]=yes, anything else=no : " -n 1 yy
	echo
	if [[ $yy =~ ^[yY]$ ]]; then
		do_it=1
	else
		echo -e "Won't proceed, I promise. Continuing on ..."
		do_it=0
	fi
}
