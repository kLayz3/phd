include q.mk
include includes.mk

# Go4 cflags needed for TFRSMap/TFOOTMap. Linking done in main Makefile tho.
CXXFLAGS += $(GO4_CXXFLAGS)
LIBS += $(GO4_LIBS)

SRC_DIR := $(shell pwd -P)
BUILD_DIR := $(SRC_DIR)/$(BUILD_EXTENSION)
DICT_DIR := $(SRC_DIR)/dict/$(BUILD)

LIBS += $(shell python3-config --embed --ldflags) \
		-lfftw3

SRC := $(wildcard $(SRC_DIR)/*.cxx)
HEADERS := $(wildcard $(SRC_DIR)/*.h)
OBJ := $(patsubst $(SRC_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(SRC))

TARGET := Structures

DICT := $(DICT_DIR)/$(TARGET)Dict.cxx
DICT_PCM := $(DICT_DIR)/$(TARGET)Dict_rdict.pcm
DICT_OBJ := $(patsubst $(DICT_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(DICT))

GEN_LIB = $(BUILD_DIR)/lib$(TARGET).so
GEN_MAP = $(BUILD_DIR)/lib$(TARGET).rootmap

# Include paths and preprocessor needed only.
RDFLAGS = $(filter -I% -D% -U%,$(CXXFLAGS)) \
		  -I../ -v3 \
		  -rml $(GEN_LIB) -rmf $(GEN_MAP)

UTIL_DIR = $(SRC_DIR)/util
UTIL_LIB = $(UTIL_DIR)/$(BUILD_EXTENSION)/libUtil.a

.PHONY: all symlink clean remove FORCE
all: $(GEN_LIB) $(DICT) symlink

# I fucked up the linking order at some point
# , so from the static I gotta load all the symbols. 🤦 */
$(GEN_LIB) : $(UTIL_LIB) $(OBJ) $(DICT_OBJ)
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)$(CXX) $(LDFLAGS) -shared -o $@ \
		-Wl,--whole-archive $< \
		-Wl,--no-whole-archive $(filter-out $<,$^) \
		$(LIBS)

$(OBJ) : $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cxx $(SRC_DIR)/%.h
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(MKDIR)
	$(Q)$(CXX) -c -o $@ $< $(CXXFLAGS)

$(DICT_OBJ) : $(DICT)
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(MKDIR)
	$(Q)$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(DICT) $(DICT_PCM) $(GEN_MAP) &: $(HEADERS) LinkDef.hh
	$(Q)$(call log,RCLING,$(notdir $(DICT)))
	@mkdir -p $(DICT_DIR) $(BUILD_DIR)
	$(Q)rootcling -f $(DICT) $(RDFLAGS) -p $^

symlink : $(GEN_LIB) $(GEN_MAP) $(DICT_PCM)
	$(Q)$(call log,SYMLINK,$(TARGET))
	$(Q)ln -sf $(DICT_PCM) $(BUILD_DIR)
	$(Q)ln -sf $(GEN_LIB) .
	$(Q)ln -sf $(GEN_MAP) .
	$(Q)ln -sf $(DICT_PCM) .

# Remove all binary/dict files for this build mode
clean: remove
	$(Q)rm -rf $(BUILD_DIR) $(DICT_DIR)
	$(MAKE) -C $(UTIL_DIR) BUILD=$(BUILD) clean

# Remove symlinks and the selected shared lib.
# Preserve objects and dict (in this build mode).
remove:
	$(Q)rm -f ./$(TARGET)Dict_rdict.pcm \
			$(BUILD_DIR)/$(TARGET)Dict_rdict.pcm \
			lib$(TARGET).so \
			lib$(TARGET).rootmap \
	    	$(GEN_LIB)
	$(MAKE) -C $(UTIL_DIR) BUILD=$(BUILD) remove

FORCE:
$(UTIL_LIB): FORCE
	$(MAKE) -C $(UTIL_DIR) BUILD=$(BUILD)
