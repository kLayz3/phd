CC=$(wildcard *.cc)
EXE =$(patsubst %.cc, %, $(CC))
AUX=$(shell root-config --cflags --evelibs --ldflags)

all: $(EXE)

% : %.cc
	g++ -o $@ $< $(AUX)
