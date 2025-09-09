OPTIMIZATION = 1

CXXFLAGS += -DPOOL_MAX_THREADS_=$(shell nproc) \
			-DPROG_PATH=\"$(shell pwd -P)\"
#CXXFLAGS += -DANALYSIS_SINGLETHREADED

ifneq ($(OPTIMIZATION),1)
CXXFLAGS += -ggdb -g3 -O0
else
CXXFLAGS += -O3 -march=native
endif

#CXXFLAGS +=  -std=c++20 
#CXXFLAGS += -Wno-cpp
