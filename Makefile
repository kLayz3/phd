CXX:=g++
SRC_DIR = src
INC_DIR = includes
BUILD_DIR = build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)
SCRIPT_DIR = scripts
GO4_SRC_DIR = $(shell pwd -P)/../go4/src
CXXFLAGS := $(shell root-config --cflags) \
	-Wall -MMD -MP -fPIC \
	-I$(INC_DIR) \
	-I$(GO4_SRC_DIR) \
	-I$(GO4SYS)/include \

include includes/common.mk

LDFLAGS := $(shell root-config --ldflags) 
LIBS := $(shell root-config --libs) \
		-L$(shell pwd -P)/includes/build \
		-Wl,-rpath,$(shell pwd -P)/includes/build \
		-lStructures \
		-L. -lGo4UserAnalysis

SRC:=$(wildcard $(SRC_DIR)/*.cc)

OBJ:=$(patsubst $(SRC_DIR)/%.cc,  $(BUILD_DIR)/%.o, $(SRC))
EXE:=$(patsubst $(SRC_DIR)/%.cc, %, $(SRC))

STRUCT_LIB = $(INC_DIR)/$(BUILD_DIR)/libStructures.so

MKDIR = mkdir -p $(@D)

.PHONY: all
all: $(EXE) $(AUX)

$(EXE): % : $(BUILD_DIR)/%.o $(STRUCT_LIB)
	ln -sf $(STRUCT_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	$(MKDIR)
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(STRUCT_LIB):
	$(MAKE) -C $(INC_DIR)

$(AUX) : $(SCRIPT_DIR)/*.bash
	./$<

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(EXE)
	$(MAKE) -C $(INC_DIR) clean
	rm -rf libStructures.so
