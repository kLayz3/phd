#!/usr/bin/env bash

help() {
	echo "This script will generate a .C macro with a corresponding entry directly as a \"void\" function."
	echo -e "Usage:\n$0 MACRO_NAME"
}

[ $# -lt 1 ] && { help; exit 1; }

name=$1
file="$name.C"

[ -f $file ] && { echo "File exists already!"; exit 2; }

echo -e "#include \"ROOT/RNTupleModel.hxx\""  >> $file 
echo -e "#include \"ROOT/RNTupleReader.hxx\"" >> $file
echo -e "#include \"ROOT/RNTupleWriter.hxx\"" >> $file
echo -e "#include \"ROOT/RCanvas.hxx\"" >> $file
echo -e "#include \"ROOT/RNTupleDS.hxx\"" >> $file
echo -e "#include \"ROOT/RDataFrame.hxx\"" >> $file
echo -e "using namespace ROOT;" >> $file
echo -e "using namespace ROOT::Experimental;\n" >> $file

echo -e "void $name(std::string fileName = \"\") {" >> $file
echo -e "\tROOT::EnableImplicitMT();" >> $file
echo -e "\n}" >> $file
nvim $file 2>/dev/null
