OPTIMIZATION = 1

ifneq ($(OPTIMIZATION),1)
CXXFLAGS += -ggdb -g3 -O0
else
CXXFLAGS += -O3 -march=native
endif

