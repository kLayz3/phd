include ../q.mk

SRC_DIR = $(shell pwd -P)
BUILD_DIR = $(SRC_DIR)/build_$(shell gcc -dumpmachine)_$(shell gcc -dumpversion)

CXXFLAGS := $(shell root-config --cflags) -Wall -MMD -MP -fPIC -I../

include ../monad/common.mk

LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --libs)

SRC := $(wildcard $(SRC_DIR)/*.cxx)
HEADERS := $(wildcard $(SRC_DIR)/*.h)

OBJ := $(patsubst $(SRC_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(SRC))

TARGET := libUtil.a
GEN_LIB = $(BUILD_DIR)/$(TARGET)

.PHONY: all symlink clean
all: $(GEN_LIB) symlink

$(GEN_LIB) : $(OBJ)
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)ar rcs $(GEN_LIB) $^

$(OBJ) : $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cxx $(SRC_DIR)/%.h
	$(Q)mkdir -p $(@D)
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(CXX) -c -o $@ $< $(CXXFLAGS)

symlink : $(GEN_LIB)
	$(Q)$(call log,SYMLINK,$(notdir $@))
	$(Q)ln -sf $(GEN_LIB) 

clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(TARGET) 
