include includes/q.mk
include includes/includes.mk

.DEFAULT_GOAL := all

BASE := $(CURDIR)
SRC_DIR = $(BASE)/src
BUILD_DIR = $(BASE)/build/$(BUILD_EXTENSION)

STRUCT_DIR = $(INC_DIR)/$(BUILD_EXTENSION)
STRUCT_LIB = $(STRUCT_DIR)/libStructures.so

CXXFLAGS += $(GO4_CXXFLAGS)
LDFLAGS += -L$(STRUCT_DIR) \
           -Wl,-rpath,$(STRUCT_DIR)
LIBS += $(GO4_LIBS)

SRC:=$(wildcard $(SRC_DIR)/*.cc)
OBJ:=$(patsubst $(SRC_DIR)/%.cc, $(BUILD_DIR)/%.o, $(SRC))
EXE:=$(patsubst $(SRC_DIR)/%.cc, $(BUILD_DIR)/%.exe, $(SRC))
LINKS:=$(patsubst $(SRC_DIR)/%.cc, %.exe, $(SRC))

.PHONY: all symlink remove clean list FORCE

all: symlink

# Always let the library's Makefile check its own prerequisites.
FORCE:

$(STRUCT_LIB): FORCE
	$(MAKE) -C $(INC_DIR) BUILD=$(BUILD)

$(EXE): $(BUILD_DIR)/%.exe : $(BUILD_DIR)/%.o $(STRUCT_LIB)
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)$(CXX) -o $@ $< $(LDFLAGS) -lStructures $(LIBS)

$(OBJ): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	@$(MKDIR)
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(CXX) -c -o $@ $< $(CXXFLAGS)

# Refresh symlinks.
symlink: $(EXE)
	$(Q)$(foreach exe, $(EXE), ln -sf $(exe) $(BASE)/$(notdir $(exe));)

# Preserve objects in this build mode.
remove:
	$(Q)rm -f $(LINKS) $(EXE)
	$(MAKE) -C $(INC_DIR) BUILD=$(BUILD) remove

clean:
	$(Q)rm -f $(LINKS)
	$(Q)rm -rf $(BUILD_DIR)
	$(MAKE) -C $(INC_DIR) BUILD=$(BUILD) clean

# h4ck1ng
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(firstword $(MAKEFILE_LIST)) : 2>/dev/null \
	  | awk -v RS= -F: '/(^|\n)# Files(\n|$$)/,/(^|\n)# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' \
	  | sort \
	  | grep -E -v -e '^[^[:alnum:]]' -e '^$@$$'
