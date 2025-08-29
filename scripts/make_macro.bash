#!/usr/bin/env bash

help() {
	echo "This script will generate a .C macro with a corresponding entry directly as a \"void\" function."
	echo -e "Usage:\n$0 MACRO_NAME"
}

[ $# -lt 1 ] && { help; exit 1; }

name=$1
file="$name.C"

[ -f $file ] && { echo "File exists already!"; exit 2; }

echo -e "void $name(const char* fileName = \"\") {\n\t\n}" >> $file
nvim $file
