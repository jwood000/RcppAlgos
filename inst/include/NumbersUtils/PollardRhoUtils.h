#pragma once

#include "PrimesPolRho.h"
#include <algorithm>
#include <vector>
#include "CppConvert/GmpxxCopy.h"

template <typename typeReturn>
void GetPrimeFactors(std::int64_t &t, std::vector<typeReturn> &factors);

/* Number of Miller-Rabin tests to run when not proving primality. */
constexpr int MR_REPS = 25;

int IsPrime(std::int64_t n);
void PollardRhoMpzT(mpz_class &n, unsigned long int a,
                    std::vector<double> &factors);
void PollardRho(std::int64_t n, std::int64_t a, std::vector<int>& factors);
