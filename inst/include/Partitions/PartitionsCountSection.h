#ifndef PARTITIONS_COUNT_SECTION_H
#define PARTITIONS_COUNT_SECTION_H

#include <cstdint>
#include <gmp.h>

// This is the maximum n such that SumSection(n) < 2^53 - 1.
// After this, we must use void SumSection(mpz_t n, mpz_t res)
constexpr int typeSwitchBnd = 328764948;

std::uint64_t SumSection(std::uint64_t n);
void SumSection(mpz_t n, mpz_t res);

#endif
