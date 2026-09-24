include includes/q.mk
include includes/includes.mk
include includes/monad/common.mk

CXX:=g++
SRC_DIR = src
BUILD_DIR = build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)
SCRIPT_DIR = scripts

CXXFLAGS += -I$(GO4SYS)/include

LIBS += -L$(shell pwd -P)/includes/build \
		-Wl,-rpath,$(INC_DIR)/build -lStructures \
		-L. -lGo4UserAnalysis

SRC:=$(wildcard $(SRC_DIR)/*.cc)

OBJ:=$(patsubst $(SRC_DIR)/%.cc,  $(BUILD_DIR)/%.o, $(SRC))
EXE:=$(patsubst $(SRC_DIR)/%.cc, %.exe, $(SRC))

STRUCT_LIB = $(INC_DIR)/$(BUILD_DIR)/libStructures.so

.PHONY: all

all: $(EXE) $(AUX)

$(EXE): %.exe : $(BUILD_DIR)/%.o $(STRUCT_LIB)
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	@$(MKDIR)
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(CXX) -c -o $@ $(CXXFLAGS) $<

$(STRUCT_LIB):
	$(MAKE) -C $(INC_DIR)

$(AUX) : $(SCRIPT_DIR)/*.bash
	$(call log,BASH,$(notdir $@))
	$(Q)./$<

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(EXE)
	$(MAKE) -C $(INC_DIR) clean
	rm -rf libStructures.so

# h4ck1ng
.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(firstword $(MAKEFILE_LIST)) : 2>/dev/null \
		| awk -v RS= -F: '/(^|\n)# Files(\n|$$)/,/(^|\n)# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' \
		| sort \
		| grep -E -v -e '^[^[:alnum:]]' -e '^$@$$'
