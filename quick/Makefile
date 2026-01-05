SRC=$(wildcard *.cc)
EXE=$(patsubst %.cc, %, $(SRC))

CFLAGS = $(shell root-config --cflags) \
	-Wall -Wextra -I../includes/

CFLAGS += -g3 -ggdb -O0 -fno-inline -fno-omit-frame-pointer

LIBS = $(shell root-config --ldflags --libs) \
	   -lROOTNTuple \
	   -L../ \
	   -Wl,-rpath,../ \
	   -lStructures

all: $(EXE)

$(EXE): % : %.cc
	g++ -o $@ $^ $(CFLAGS) $(LIBS)
