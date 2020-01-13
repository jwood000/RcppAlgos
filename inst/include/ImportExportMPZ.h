#ifndef IMPORT_EXPORT_MPZ_H
#define IMPORT_EXPORT_MPZ_H

#include <gmp.h>
#include <Rcpp.h>

// Functions for importing/exporting and converting SEXPs to mpz_t

void createMPZArray(SEXP input, mpz_t *myVec, std::size_t vecSize, 
                     const std::string &nameOfObject, bool negPoss = false);

int myRaw(char* raw, mpz_t value, std::size_t totals);

#endif
