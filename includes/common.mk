OPTIMIZATION = 1

CXXFLAGS += -DPOOL_MAX_THREADS_=$(shell nproc) \
			-DPROG_PATH=\"$(shell pwd -P)\"

ifneq ($(OPTIMIZATION), 1)

CXXFLAGS += -ggdb3 -g3 -O0 -fno-omit-frame-pointer -fno-inline
export ASAN_OPTIONS=detect_leaks=1,strict_string_checks=1,alloc_dealloc_mismatch=1
export MALLOC_CHECK_=3
export GLIBC_TUNABLES=glibc.malloc.tcache_count=0

else

CXXFLAGS += -O3 -march=native

endif
