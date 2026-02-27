#!/usr/bin/env bash

script_dir=$(realpath $(dirname -- $0)

help() {
	echo "This script will generate a .cc source template."
	echo -e "Usage:\n$0 STEP MACRONAME"
	echo -e "Where STEP is one of: map cal hit"
}

[ $# -lt 2 ] && { help; exit 2; }

step=$1
[[ $step != "map" && $step != "cal" && $step != "hit" ]] && { echo -e "Step input wrong\n"; help; exit 3; }

name=$2
[[ ! $name =~ ^[[:alpha:]].*[[:alnum:]]$ ]] && { echo -e "Script name must start with letter [a-zA-Z] and end with alphanumeric\n"; help; exit 4; }
file=$step/$name.cc

mkdir -p $step
[ -f $file ] && { echo "File exists already!"; exit 2; }

echo -e "#include \"../libs.hh\""  > $file 

for x in $(find ../includes -maxdepth 1 -type f -name "*.h" -o -name "*.hxx"); do
	echo -e "#include \"$x\"" >> $file
done

echo -e "\nnamespace $name {\n" >> $file
echo -e "extern const char* help;\n"
echo -e "void $name(int argc, char** argv) {" >> $file
echo -e "using namespace indicators;" >> $file
echo -e "auto& def_msg = CMDLineParser::Mandatory::DefMessage";
echo -e "CMDLineParser::Mandatory::SetDefMessage(help)";
echo -e "\tROOT::EnableImplicitMT();" >> $file
echo -e "\tTApplication* app = new TApplication(\"myApp\", 0, 0);" >> $file
echo -e "\n\tprintf(\"End of script: $name.\");" >> $file
echo -e "\tapp->Run();" >> $file
echo -e "\n}\n" >> $file
echo -e "const char* help = \n\"Usage: ./$name ARGS\"" >> $file
echo -e "}" >> $file # end of namespace


nvim $file 2>/dev/null
