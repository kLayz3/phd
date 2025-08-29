#pragma once

#include <chrono>
#include <cstdio>
#include "TTree.h"
#include "TGraph.h"
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>
#include "libs.hh"
#include "indicators.hh"

inline void PrintProgress(indicators::ProgressBar& bar,	u64 n_entry, u64 max_entries, u64 step = 100) noexcept {
	static u64 n_entry_called = 0;
	if(n_entry - n_entry_called < step) return;
	bar.set_progress( (n_entry*100) / max_entries );

	n_entry_called = n_entry;
}

/* ------------------------- */
inline uint64_t SortEntries(uint64_t& firstEvent, uint64_t& maxEvents, TTree* h101) {
	firstEvent = std::min(firstEvent, (uint64_t)h101->GetEntries());
    uint64_t n = (maxEvents==0 || firstEvent+maxEvents > (uint64_t)h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
    maxEvents = n - firstEvent;
	return n;
}

/* ------------------------- */
template<class A, class B>
inline void TGraphFromVector(TGraph& t, std::vector<A> vX, std::vector<B> vY) {
	auto s = std::min(vX.size(), vY.size());
	for(auto i=0; i<s; ++i) {
		t.AddPoint((double)vX[i], (double)vY[i]);
	}
}

/* ------------------------- */
template<class T>  
inline void QuickSwap(std::vector<T>& v, int i, int j) {
	if(!v.size()) return;
    std::swap(v[i], v[j]);
}

/* ------------------------- */
template<class T>  
inline void QuickErase(std::vector<T> &v, int i) {
    std::swap(v[i], v.back());
    v.pop_back();
}

/* ------------------------- */
template<class T>
inline void ReleaseMalloc(T x) {
	free(x);
}
template<class T, class... Args>
inline void ReleaseMalloc(T x, Args... args) {
	free(x); ReleaseMalloc(args...);
}

/* ------------------------- */
template<class T>
inline T rround(double x) noexcept { return static_cast<T>(x + 0.5); }

using std::chrono::duration_cast;
using std::chrono::seconds;
using std::chrono::milliseconds;

#define timeNow() std::chrono::high_resolution_clock::now()

struct TimePoint {
	std::chrono::high_resolution_clock::time_point t;
	std::string tag;
	TimePoint() : t(timeNow()), tag("?") {};
	TimePoint(std::string s) : t(timeNow()), tag(s) {};
};

enum TimingVariant {
	kMINUTE,
	kSECOND,
	kMILLISECOND,
	kMICROSECOND,
};

template<TimingVariant E = kMILLISECOND>
void PrintElapsed(const TimePoint& end, const TimePoint& start) {
	using ms = std::chrono::milliseconds;
	using us = std::chrono::microseconds;
	
	double elapsed; const char* mode;
	if constexpr(E == kMINUTE) {
		elapsed = static_cast<double>(std::chrono::duration_cast<ms>(end.t - start.t).count()) / 60'000;
		mode = "min";
	}
	else if constexpr(E == kSECOND) {
		elapsed = static_cast<double>(std::chrono::duration_cast<ms>(end.t - start.t).count()) / 1'000;
		mode = "sec";
	}
	else if constexpr(E == kMILLISECOND) {
		elapsed = static_cast<double>(std::chrono::duration_cast<us>(end.t - start.t).count()) / 1'000;
		mode = "ms";
	}
	else if constexpr(E == kMICROSECOND) {
		elapsed = static_cast<double>(std::chrono::duration_cast<us>(end.t - start.t).count());
		mode = "us";
	}
    printf("Elapsed time from " EMPH1(%s) " to " EMPH1(%s) ": " EMPH(%.1f) " %s\n", start.tag.c_str(), end.tag.c_str(), elapsed, mode);
}

template<TimingVariant E = kMILLISECOND>
inline void PrintElapsed(const std::vector<TimePoint>& v) {
	if(v.size() < 2) return;
	PrintElapsed<E>(v[v.size()-1], v[v.size() - 2]);
}


inline void __concat_impl__(std::ostringstream& ) {}

template<typename T, typename... Args>
void __concat_impl__(std::ostringstream& oss, T&& first, Args&&... args) {
	oss << std::forward<T>(first);
	__concat_impl__(oss, std::forward<Args>(args)... );
}

/**
 * Concatenates bunch of arguments which can be lvalues, statics, etc. and returns an owned std::string. 
 */
template<typename... Args>
std::string sstrcat(Args&&... args) {
	std::ostringstream oss;
	__concat_impl__(oss, std::forward<Args>(args)... );
	return oss.str();
}

namespace util {
	template<typename T>
	struct is_an_array : std::false_type {};

	template<typename T, size_t N>
	struct is_an_array<T[N]> : std::true_type { using value_type = T; static constexpr size_t size = N; };

	template<typename T, size_t N>
	struct is_an_array<std::array<T,N>> : std::true_type { using value_type = T; static constexpr size_t size = N; };

	template<typename T>
	inline constexpr bool is_an_array_v = is_an_array<T>::value;
}

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>,
	typename std::enable_if<util::is_an_array_v<U>>::type* = nullptr>
auto median(const T& arr) {
	constexpr std::size_t N = util::is_an_array<U>::size;
	if constexpr(N % 2)
		return arr[(N-1) / 2];
	else 
		return ( arr[(N-1) / 2] + arr[N/2] ) / 2;
}

template<typename T,
	typename U = std::remove_cv_t<std::remove_reference_t<T>>,
	typename std::enable_if<util::is_an_array_v<U>>::type* = nullptr>
constexpr int FindIndex(const T& arr, const typename util::is_an_array<U>::value_type& val) {
	constexpr std::size_t N = util::is_an_array<U>::size;
	for(int i=0; i < (int)N; ++i)
		if(arr[i] == val) return i;
	return -1;
}
