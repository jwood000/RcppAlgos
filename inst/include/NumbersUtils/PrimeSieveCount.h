#ifndef PRIME_SIEVE_COUNT_H
#define PRIME_SIEVE_COUNT_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP PrimeSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads,
                       SEXP RmaxCores, SEXP RmaxThreads);
    SEXP PrimeCountCpp(SEXP Rn, SEXP RNumThreads,
                       SEXP RmaxThreads);
}

#endif
