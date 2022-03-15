#ifndef MOTLEY_PRIMES_H
#define MOTLEY_PRIMES_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP MotleyContainer(SEXP Rb1, SEXP Rb2, SEXP isEuler, SEXP RNamed,
                         SEXP RNumThreads, SEXP maxThreads);
}

#endif
