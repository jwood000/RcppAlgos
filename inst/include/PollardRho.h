#ifndef POLLARDRHO_R
#define POLLARDRHO_R 1

#include <stdint.h>

/**
 * Get Prime Factorization
 * t: number to factorize
 * factors [out]: the list of factors
 *
 * Note: this is adapted from demo "factorize.c" file from gmplib
 */
void getPrimefactors (int64_t& t, std::vector<int64_t>&  factors);

#endif
