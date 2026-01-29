#!/bin/bash

# This scripts will generate .clangd file in the main folder, to allow proper LSP linting of
# ROOT and Go4 libraries.

script_dir=$(dirname -- $(readlink -f -- $0))
CLANGD=$script_dir/../.clangd
INC_DIR=$script_dir/../includes

echo \
"CompileFlags:
  Add:
    - -I$(root-config --incdir)
    - -I$INC_DIR
    - -std=c++17
  Remove: [-std=*]
" > $CLANGD

#CERES_DIR=ceres-solver
	#- -I$INC_DIR/$CERES_DIR/include 
	#- -I$INC_DIR/$CERES_DIR/third_party/abseil-cpp 
	#- -I$INC_DIR/$CERES_DIR/build/include
#Compiler: /usr/bin/g++
