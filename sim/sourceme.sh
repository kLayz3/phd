. /etc/os-release

CPP_VER=17
FAIRSOFT_BASE_DIR=/cvmfs/fairsoft.gsi.de/debian${VERSION_ID}/fairsoft/cpp$CPP_VER

# Just take latest.
FAIRSOFT_DIR=$(ls -t1 $FAIRSOFT_BASE_DIR | head -1)

# Check if `thisroot.sh` exists. It should
THISROOT="$FAIRSOFT_BASE_DIR/$FAIRSOFT_DIR/bin/thisroot.sh"
if [ ! -f $THISROOT ]; then
    echo -e "\nIdentified host: $(hostname), but attempting to source fairsoft root failing?"
    echo "Source file not found: $THISROOT"
    exit 2;
fi

source $THISROOT
