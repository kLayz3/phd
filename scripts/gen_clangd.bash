#!/bin/bash

# This scripts will generate .clangd file in the main folder, to allow proper LSP linting of
# ROOT and Go4 libraries.

script_dir=$(dirname -- $(readlink -f -- $0))
CLANGD=$script_dir/../.clangd
INC_DIR=$script_dir/../includes

rm -f $CLANGD && touch $CLANGD

echo \
"CompileFlags:
	Add:
		- \"-I$(root-config --incdir)\"
		- \"-I$INC_DIR\"
		- \"-I$script_dir/../\"
" >> $CLANGD
