#ifndef POLLARD_RHO_H
#define POLLARD_RHO_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP PollardRhoContainer(SEXP Rv, SEXP RNamed,
                             SEXP RbPrimeFacs, SEXP RbAllFacs,
                             SEXP RNumThreads, SEXP RmaxThreads);
}

#endif
