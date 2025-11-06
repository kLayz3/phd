GO4_SRC_PATH = $(shell pwd -P)/../../go4/src

CFLAGS = $(shell root-config --cflags) -I../includes \
	   -I$(GO4SYS)/include \
	   -I$(GO4_SRC_PATH) \
	   -g -ggdb \
	   -DGO4_SRC_PATH=\"$(GO4_SRC_PATH)\"

LIBS = $(shell root-config --libs --ldflags) \
	 -lROOTNTuple -lROOTNTupleUtil

LIBS += -L$(GO4_SRC_PATH) \
		-Wl,-rpath,$(GO4_SRC_PATH) -lGo4UserAnalysis \

TAR = gen_rntuple gen_rntuple_frs

all: $(TAR)

$(TAR): % : %.cxx
	g++ -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(TAR)
