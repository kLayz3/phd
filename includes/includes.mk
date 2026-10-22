# The realpath of this specific file.
INC_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

include $(INC_DIR)/monad/common.mk

# Build directory always comes suffixed with the tag from monad's common.
BUILD_EXTENSION := build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)/$(BUILD)

# Flags..
CXXFLAGS += $(shell root-config --cflags) -Wall -fPIC \
           -I$(INC_DIR)

LDFLAGS += $(shell root-config --ldflags)
LIBS += $(shell root-config --libs) -lROOTNTuple \
        -Wl,-rpath,'$$ORIGIN'

MKDIR = mkdir -p $(@D)

CXX = g++

# Go4 related stuff.

GO4_SRC_DIR := $(INC_DIR)/../../go4/src

GO4_CXXFLAGS := -I$(GO4_SRC_DIR) \
                -I$(GO4SYS)/include

GO4_LIBS := -Wl,-rpath,$(GO4_SRC_DIR) \
            -L$(GO4_SRC_DIR) \
            -lGo4UserAnalysis
