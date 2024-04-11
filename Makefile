GCC:=g++
SRC_DIR:=src
INC_DIR:=includes
BUILD_DIR:=build

CFLAGS:=-c -g -Wall
CFLAGS+=$(shell root-config --cflags)\
		$(shell root-config --auxcflags)\
		-I$(INC_DIR) -I.

LDFLAGS:=$(shell root-config --ldflags)
ROOTLIBS:=$(shell root-config --libs)

SRC:=$(wildcard $(SRC_DIR)/*.cc)
INC:=$(wildcard $(INC_DIR)/*.cxx)

OBJ:=$(patsubst $(SRC_DIR)/%.cc,  $(BUILD_DIR)/%.o,   $(SRC))\
	 $(patsubst $(INC_DIR)/%.cxx, $(BUILD_DIR)/%.oxx, $(INC))

EXE:=$(patsubst $(SRC_DIR)/%.cc, %, $(SRC))

MKDIR = mkdir -p $(@D)

.PHONY: all
all: $(EXE)

$(EXE) : $(OBJ)
	$(GCC) $(LDFLAGS) $(BUILD_DIR)/$@.o $(BUILD_DIR)/*.oxx -o $@ $(ROOTLIBS)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cc
	$(MKDIR)
	$(GCC) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.oxx : $(INC_DIR)/%.cxx
	$(MKDIR)
	$(GCC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) *~
	rm -f $(EXE)
