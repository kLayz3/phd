#pragma once

#define FOOT25 25
#define FOOT23 23
#define FOOT22 22
#define FOOT21 21
#define FOOT20 20
#define FOOT19 19
#define FOOT17 17
#define FOOT10 10

#include <unordered_set>
#include <unordered_map>
/* Table of bad strips, per detector (searched by the same index as in the ROOT file) */
typedef std::unordered_map<int, std::unordered_set<int>> BadStripTable;

//extern BadStripTable bad_strips = {
//};
