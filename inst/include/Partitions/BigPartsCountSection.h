#ifndef BIG_PARTS_COUNT_SECTION_H
#define BIG_PARTS_COUNT_SECTION_H

#include <gmp.h>

// This is the maximum n such that SumSection(n) < 2^53 - 1.
// After this, we must use void SumSection(mpz_t n, mpz_t res)
constexpr int typeSwitchBnd = 328764948;

void SumSection(mpz_t n, mpz_t res);

#endif
