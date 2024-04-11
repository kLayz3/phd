#pragma once

#include <cstdio>
#include "TTree.h"
#include "TGraph.h"
#include <vector>

#define PBARSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBARW 60

inline void PrintProgress(float percentage) {
    int val = (int)(percentage*100);
    int lpad =(int)(percentage*PBARW);
    int rpad = PBARW-lpad;
    printf("\r%3d%% [%.*s%*s]",val,lpad,PBARSTR,rpad,"");
    fflush(stdout);
}

/* ------------------------- */
inline uint64_t SortEntries(uint64_t& firstEvent, uint64_t& maxEvents, TTree* h101) {
	firstEvent = std::min(firstEvent, (uint64_t)h101->GetEntries());
    uint64_t n = (maxEvents==0 || firstEvent+maxEvents > h101->GetEntries()) ? (h101->GetEntries()) : (firstEvent+maxEvents);
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
