#ifndef DIV_NUM_SIEVE_H
#define DIV_NUM_SIEVE_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP DivNumSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RbDivSieve,
                        SEXP RisNamed, SEXP RNumThreads,
                        SEXP RmaxThreads);
}

#endif
