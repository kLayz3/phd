BASE := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

include $(BASE)/../includes/q.mk
include $(BASE)/../includes/includes.mk

# Keep Go4 include paths out of these programs.
# Basically, adding this flag to the Go4 includes freaks the compiler out, lol.
CXXFLAGS += -Wextra -I$(BASE)

BUILD_DIR := $(BASE)/$(BUILD_EXTENSION)
STRUCT_DIR := $(INC_DIR)/$(BUILD_EXTENSION)
STRUCT_LIB := $(STRUCT_DIR)/libStructures.so

SUBDIRS := map cal hit
SRC := $(BASE)/run.cc
SRC += $(foreach subdir,$(SUBDIRS), $(wildcard $(BASE)/$(subdir)/*.cc))

EXE := $(patsubst $(BASE)/%.cc,$(BUILD_DIR)/%,$(SRC))
LINKS := $(SRC:.cc=)

LDFLAGS += -L$(STRUCT_DIR) -Wl,-rpath,$(STRUCT_DIR)

.PHONY: all symlink remove clean
all: symlink

# Build the selected library if it is missing.
$(STRUCT_LIB):
	$(MAKE) -C $(INC_DIR) BUILD=$(BUILD)

$(EXE): $(BUILD_DIR)/%: $(BASE)/%.cc $(STRUCT_LIB)
	$(Q)$(MKDIR)
	$(Q)$(call log,CXX,$(patsubst $(BUILD_DIR)/%,%,$@))
	$(Q)$(CXX) -o $@ $< $(CXXFLAGS) \
	    $(LDFLAGS) -lStructures $(LIBS)

# Retain the existing progs: run, map/name, cal/name, hit/name.
symlink: $(EXE)
	$(Q)$(foreach exe,$(EXE),ln -sf $(exe) $(BASE)/$(patsubst $(BUILD_DIR)/%,%,$(exe));)

# Keep objects for relinking...
remove:
	$(Q)rm -f $(LINKS)

clean:
	$(Q)rm -f $(LINKS)
	$(Q)rm -rf $(BUILD_DIR)

# Clean never cleans or alters the main files from $(BASE)../
