include q.mk

SRC_DIR = $(shell pwd -P)
BUILD_DIR = $(SRC_DIR)/build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)
DICT_DIR = $(SRC_DIR)/dict

GO4_SRC_DIR = $(shell pwd -P)/../../go4/src
CXXFLAGS := $(shell root-config --cflags) -Wall -MMD -MP -fPIC \
	-I$(SRC_DIR)/core \
	-I$(GO4_SRC_DIR) \
	-I$(GO4SYS)/include

include monad/common.mk

LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --libs) -lROOTNTuple \
		-Wl,-rpath,'$$ORIGIN'

# Normally on lustre.hpc the libfftw3 is fucked so I found
# one normally on a /usr/lib that works and ported it over to $HOME/.local/lib
# This is a very handwavy hack.
HOSTNAME = $(shell hostname)
ifneq ($(findstring lustre,$(HOSTNAME)),)
	LIBS += -Wl,--no-as-needed \
		-L$(HOME)/.local/lib \
		-Wl,-rpath,$(HOME)/.local/lib \
		-lfftw3 \
		-Wl,--as-needed
else
	LIBS += -lfftw3
endif

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
		  -I$(GO4_SRC_DIR) \
		  -Wall -MMD -MP -rml $(GEN_LIB) -rmf $(GEN_MAP)

UTIL_DIR = $(SRC_DIR)/util
UTIL_LIB = $(UTIL_DIR)/libUtil.a

.PHONY: all symlink clean
all: $(GEN_LIB) $(DICT) symlink

# I fucked up the linking order at some point
# , so from the static I gotta load all the symbols. 🤦 */ 
$(GEN_LIB) : $(UTIL_LIB) $(OBJ) $(DICT_OBJ) 
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)$(CXX) $(LDFLAGS) -shared -o $@ \
		-Wl,--whole-archive $< \
		-Wl,--no-whole-archive $(filter-out $<,$^) \
		$(LIBS)	

$(OBJ) : $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cxx $(SRC_DIR)/%.h $(DICT_OBJ)
	$(Q)$(call log,CXX,$(notdir $@))
	@mkdir -p $(@D)
	$(Q)$(CXX) -c -o $@ $< $(CXXFLAGS)

$(DICT_OBJ) : $(DICT)
	$(Q)$(call log,CXX,$(notdir $@))
	@mkdir -p $(@D)
	$(Q)$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(DICT) : $(HEADERS) LinkDef.hh 
	$(Q)$(call log,RCLING,$(notdir $@))
	@mkdir -p $(DICT_DIR) $(BUILD_DIR)
	$(Q)rootcling -f $@ $(RDFLAGS) -p $^

$(UTIL_LIB):
	$(MAKE) -C $(UTIL_DIR)

symlink : $(DICT) $(GEN_LIB)
	$(Q)ln -sf $(DICT_PCM) $(BUILD_DIR)
	$(Q)ln -sf $(GEN_LIB)  ../
	$(Q)ln -sf $(GEN_MAP)  ../
	$(Q)ln -sf $(DICT_PCM) ../

clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(DICT_DIR)
	$(MAKE) -C $(UTIL_DIR) clean
