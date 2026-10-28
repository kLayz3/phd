# The realpath of this specific file.
INC_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

include $(INC_DIR)/monad/common.mk

CXX = g++

COMPILER_TAR := $(shell $(CXX) -dumpmachine)
COMPILER_VER := $(shell $(CXX) -dumpfullversion)

# Same hosts of identical GCC version can compile down to same extension,
# even tho their arch'es and compiler tags differ. Get them here, hopefully, inside $NATIVE_ARCH var.
# The text-replacement part comes due to the stupid fact how GCC formats the output, e.g.:
# g++ -march=native -Q --help=target 2>/dev/null | sha256sum
# The output becomes different based on the terminal size... Hehehe.
# Also then, do not care about helper text and only take lines which are:
# -feature value     and hash that.
ifndef NATIVE_ARCH
NATIVE_ARCH := $(shell opts=$$(LC_ALL=C COLUMNS=300 \
			   $(CXX) -march=native -Q --help=target 2>/dev/null) && \
			   opts=$$(printf '%s' "$$opts" \
			   		   | perl -ne 'print "$$1 $$2\n" if /^\s*(-m\S+)\s+(\S+)\s*$$/') \
					&& test -n "$$opts" \
					&& { printf '%s\n' "$$opts" | LC_ALL=C sort | sha256sum | cut -c1-15; } \
					|| { echo "NULL"; })
endif

ifeq ($(NATIVE_ARCH),NULL)
  $(error Could not determine native CPU configuration)
endif

export NATIVE_ARCH

# Build directory always comes suffixed with the $(BUILD) tag from monad's common.
BUILD_EXTENSION := build_$(COMPILER_TAR)_$(COMPILER_VER)_$(NATIVE_ARCH)/$(BUILD)

# ROOT system tag,.. stupid edgecases with . and slash. LOL.
ROOT_TAG := $(subst /,-,$(subst .,-,$(shell root-config --version)))

# Flags..
CXXFLAGS += $(shell root-config --cflags) -Wall -fPIC \
           -I$(INC_DIR)

LDFLAGS += $(shell root-config --ldflags)
LIBS += $(shell root-config --libs) -lROOTNTuple \
        -Wl,-rpath,'$$ORIGIN'

MKDIR = mkdir -p $(@D)

# Go4 related stuff.

GO4_SRC_DIR := $(INC_DIR)/../../go4/src

GO4_CXXFLAGS := -I$(GO4_SRC_DIR) \
                -I$(GO4SYS)/include

GO4_LIBS := -Wl,-rpath,$(GO4_SRC_DIR) \
            -L$(GO4_SRC_DIR) \
            -lGo4UserAnalysis
