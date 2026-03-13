#pragma once

#include "PolyFitter.hxx"

/* This hack is used to instantiate the template for small-ish N's.
 * This has the advantage that then these procedures can be used in the ROOT macros.
 *
 * Problem is that Eigen is heavily optimized from rank 2 onward for algorithms such
 * as colPivHouseholderQr() decomposition and CLING when it compiles the source of macro
 * won't enable Eigen optimizations, and you will catch the weirdest segfault in your life.
 *
 * Don't ask me why, but I was swearing in 4 different languages debugging this... */
