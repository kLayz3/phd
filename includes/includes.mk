INC_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
GO4_SRC_DIR := $(INC_DIR)../../go4/src

CXXFLAGS := $(shell root-config --cflags) -Wall -MMD -MP -fPIC \
	-I$(INC_DIR) \
	-I$(GO4_SRC_DIR) \
	-I$(GO4SYS)/include

LDFLAGS := $(shell root-config --ldflags)

LIBS := $(shell root-config --libs) -lROOTNTuple \
		-Wl,-rpath,'$$ORIGIN' \

MKDIR = mkdir -p $(@D)
