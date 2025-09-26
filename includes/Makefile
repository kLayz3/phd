SRC_DIR = $(shell pwd -P)
BUILD_DIR = $(SRC_DIR)/build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)
DICT_DIR = $(SRC_DIR)/dict

CXXFLAGS := $(shell root-config --cflags) -Wall -MMD -MP -fPIC

include common.mk

LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --libs) -lROOTNTuple

SRC := $(wildcard $(SRC_DIR)/*.cxx)
HEADERS := $(wildcard $(SRC_DIR)/*.h)

OBJ := $(patsubst $(SRC_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(SRC))

TARGET := Structures
DICT := $(DICT_DIR)/$(TARGET)Dict.cxx
DICT_PCM := $(DICT_DIR)/$(TARGET)Dict_rdict.pcm

DICT_OBJ := $(patsubst $(DICT_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(DICT))
GEN_LIB = $(BUILD_DIR)/lib$(TARGET).so
GEN_MAP = $(BUILD_DIR)/lib$(TARGET).rootmap

RDFLAGS = -I./ -I../ -v3 \
		  -Wall -MMD -MP -rml $(GEN_LIB) -rmf $(GEN_MAP)

.PHONY: all symlink clean
all: $(GEN_LIB) $(DICT) symlink

$(GEN_LIB) : $(OBJ) $(DICT_OBJ)
	$(CXX) $(LDFLAGS) -shared -o $@ $^ $(LIBS)	

$(OBJ) : $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cxx $(SRC_DIR)/%.h $(DICT_OBJ)
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(DICT_OBJ) : $(DICT)
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(DICT) : $(HEADERS) LinkDef.hh 
	@mkdir -p $(DICT_DIR) $(BUILD_DIR)
	rootcling -f $@ $(RDFLAGS) -p $^

symlink : $(DICT) $(GEN_LIB)
	ln -sf $(DICT_PCM) $(BUILD_DIR)
	ln -sf $(GEN_LIB) ../
	ln -sf $(GEN_MAP) ../

clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(DICT_DIR)
