BUILD_DIR=build
SRC=$(wildcard *.cc)
EXE=$(patsubst %.cc, $(BUILD_DIR)/%, $(SRC))

CXXFLAGS = $(shell root-config --cflags) \
	-Wall -Wextra -I../includes/

include ../includes/monad/common.mk

LIBS = $(shell root-config --ldflags --libs) \
	   -lROOTNTuple \
	   -L../ \
	   -Wl,-rpath,../ \
	   -lStructures

all: $(EXE)

$(EXE): $(BUILD_DIR)/% : %.cc
	mkdir -p $(@D)
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)
