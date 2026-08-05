BASE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
BASE := $(patsubst %/,%,$(BASE))

include $(BASE)/../includes/monad/common.mk
include $(BASE)/../includes/q.mk

CXXFLAGS += $(shell root-config --cflags) \
	-Wall -Wextra -I$(MONAD_DIR)../ -I./ \

LIBS += $(shell root-config --libs)
LDFLAGS += $(shell root-config --ldflags) 

STEP := map cal hit

SRC := $(BASE)/run.cc
SRC += $(wildcard $(addsuffix /*.cc, $(STEP)))

EXE := $(patsubst %.cc, %, $(SRC))

LIBS += -lROOTNTuple \
		-L../ \
		-Wl,-rpath,../ \
		-lStructures

.PHONY: all clean
all: $(EXE)

% : %.cc
	$(Q)$(call log,CXX,$(notdir $@))
	$(Q)g++ -o $@ $< $(CXXFLAGS) $(LDFLAGS) $(LIBS)

clean:
	$(Q)$(call log,RM,$(notdir $@))
	$(Q)rm -rf $(EXE) 2>/dev/null

