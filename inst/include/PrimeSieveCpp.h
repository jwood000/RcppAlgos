#ifndef PRIME_SIEVE_H
#define PRIME_SIEVE_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP PrimeSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads,
                       SEXP RmaxCores, SEXP RmaxThreads);
}

#endif
