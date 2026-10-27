include ../q.mk
include ../includes.mk

.DEFAULT_GOAL := all

SRC_DIR = $(shell pwd -P)
BUILD_DIR = $(SRC_DIR)/$(BUILD_EXTENSION)

CXXFLAGS += -I$(INC_DIR)/pybind11/include \
            $(shell python3-config --includes)

SRC := $(wildcard $(SRC_DIR)/*.cxx)
HEADERS := $(wildcard $(SRC_DIR)/*.h)
OBJ := $(patsubst $(SRC_DIR)/%.cxx, $(BUILD_DIR)/%.o, $(SRC))

TARGET := libUtil.a
GEN_LIB = $(BUILD_DIR)/$(TARGET)

.PHONY: all symlink clean remove
all: $(GEN_LIB) symlink

$(GEN_LIB) : $(OBJ)
	$(Q)$(call log,LINK,$(notdir $@))
	$(Q)ar rcs $(GEN_LIB) $^

$(OBJ) : $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cxx $(SRC_DIR)/%.h
	$(Q)$(MKDIR)
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)$(CXX) -c -o $@ $< $(CXXFLAGS)

symlink : $(GEN_LIB)
	$(Q)$(call log,SYMLINK,$(notdir $@))
	$(Q)ln -sf $(GEN_LIB)

clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(TARGET)

remove:
	rm -rf $(TARGET)
