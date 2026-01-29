#!/usr/bin/env bash

help() {
	echo "This script will generate a .cc source template."
	echo -e "Usage:\n$0 MACRO_NAME"
}

[ $# -lt 1 ] && { help; exit 1; }

name=$1
file="$name.cc"

[ -f $file ] && { echo "File exists already!"; exit 2; }

echo -e "#include \"libs.hh\""  >> $file 

for x in $(find ../includes -maxdepth 1 -type f -name "*.h" -o -name "*.hxx"); do
	echo -e "#include \"$x\"" >> $file
done

echo -e "\nusing namespace ROOT;" >> $file
echo -e "using namespace ROOT::Experimental;\n" >> $file
echo -e "using namespace std;"

echo -e "int main(int argc, char* argv[]) {" >> $file
echo -e "\tROOT::EnableImplicitMT();" >> $file
echo -e "\tTApplication* app = new TApplication(\"myApp\", 0, 0);" >> $file
echo -e "\n\tprintf(\"End of main.\");" >> $file
echo -e "\tapp->Run();" >> $file
echo -e "\treturn 0;\n}" >> $file

nvim $file 2>/dev/null
