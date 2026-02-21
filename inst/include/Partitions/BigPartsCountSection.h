#pragma once

#include <gmpxx.h>
#include <vector>

// This is the maximum n such that SumSection(n) < 2^53 - 1.
// After this, we must use void SumSection(mpz_class n, mpz_class res)
constexpr int typeSwitchBnd = 328764948;

void SumSection(const mpz_class &n, mpz_class &res);
void ResetP2D(std::vector<std::vector<mpz_class>> &p2d);
