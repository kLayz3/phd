include common.mk

SRC := $(BASE)/main.cc
OBJ := $(patsubst $(BASE)/%.cc, $(BASE)/$(BUILD)/%.o, $(SRC_N))
EXE := $(BASE)/run

STEP := map cal hit

LIBS = $(shell root-config --ldflags --libs) \
	   -lROOTNTuple \
	   -L../ \
	   -Wl,-rpath,../ \
	   -lStructures

ALIB = $(patsubst %, $(BASE)/%/$(BUILD)/%.a, $(STEP))

.PHONY: all
all: $(EXE)

$(EXE): $(OBJ) $(ALIB)
	$(Q)echo "  LD $(notdir $@) "
	$(Q)g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

$(OBJ): $(BASE)/$(BUILD)/%.o : $(BASE)/%.cc
	$(Q)echo "  CXX $(notdir $@)  "
	$(Q)g++ -c -o $@ $< $(CXXFLAGS)

$(ALIB):
	$(Q)echo "  MAKE $(notdir $@)  "
	$(Q)$(MAKE) -C $(dir $@)..

.PHONY: clean
clean:
	$(Q)rm -rf $(BASE)/$(BUILD) $(BASE)/$(EXE) 2>/dev/null

.PHONY: clean-all $(addprefix clean-, $(STEP))

clean-all: $(addprefix clean-,$(STEP))

clean-%:
	$(Q)echo "  MAKE CLEAN $*  "
	$(Q)$(MAKE) -C $(BASE)/$* clean

