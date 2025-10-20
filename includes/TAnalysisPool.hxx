#pragma once

#include "libs.hh"
#include "TAnalysisProcess.hxx"

/* Passing in this define to purposefully undef
* the singlethreaded-ness. */
#if defined(ANALYSIS_MULTITHREADED)
#	undef ANALYSIS_SINGLETHREADED
#endif

/* Passed from the build tool. */
#if !defined(ANALYSIS_SINGLETHREADED)
#	if !defined(POOL_MAX_THREADS_)
#		define POOL_MAX_THREADS_ 10
#	endif
#else
#	warning "Running single-threaded.."
#	define POOL_MAX_THREADS 1
#endif

template <
	u32 N, u32 NSlice,
	typename... Ts
> struct TAnalysisPool final {
	static_assert(N <= POOL_MAX_THREADS_, 
		"TAnalysisProcess template instantiated with over-the-top capacity: " _TO_STRING(POOL_MAX_THREADS_) );
	static_assert(N != 0, "TAnalysisProcess template size == 0?");

private:
	std::array <
		TAnalysisProcess<Ts...>, N
	> pool;
};
