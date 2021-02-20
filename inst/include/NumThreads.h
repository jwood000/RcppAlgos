#ifndef NUM_THREADS_H
#define NUM_THREADS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP cpp11GetNumThreads();
}

#endif
