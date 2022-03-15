#ifndef CARTESIAN_CONTAINER_H
#define CARTESIAN_CONTAINER_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP ComboGridCpp(SEXP RList, SEXP RIsRep);
}

#endif
