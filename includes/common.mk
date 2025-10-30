OPTIMIZATION = 0

CXXFLAGS += -DPOOL_MAX_THREADS_=$(shell nproc) \
			-DPROG_PATH=\"$(shell pwd -P)\"

ifneq ($(OPTIMIZATION), 1)
CXXFLAGS += -ggdb -g3 -O0 -DANALYSIS_SINGLETHREADED
else
CXXFLAGS += -O3 -march=native
endif

#CXXFLAGS += -g -fno-omit-frame-pointer -fno-inline # For callgrind.
#CXXFLAGS +=  -std=c++20 
#CXXFLAGS += -Wno-cpp
