CXX ?= g++
ROOT_CONFIG ?= root-config

ROOT_CXXFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs --ldflags)
CXXFLAGS ?= -O2 -g
LDLIBS += -lROOTNTuple -lROOTNTupleUtil

TAR = gen_rntuple

.PHONY: all clean
all: $(TAR)

$(TAR): % : %.cxx
	$(CXX) $(CPPFLAGS) $(ROOT_CXXFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(ROOT_LIBS) $(LDLIBS)

clean:
	rm -f $(TAR)
