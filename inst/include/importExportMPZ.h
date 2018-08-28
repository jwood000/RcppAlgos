#ifndef IMPORT_EXPORT_MPZ_GMP_H
#define IMPORT_EXPORT_MPZ_GMP_H

#include <gmp.h>
#include <Rcpp.h>

/**
 * Functions for importing/exporting and converting SEXPs to mpz_t
 * 
 */
void createMPZArray (SEXP v, mpz_t myVec[], unsigned long int sizevec);

int myRaw (char* raw, mpz_t value, unsigned long int totals);

#endif
