#ifndef POLLARD_RHO_UTILS_H
#define POLLARD_RHO_UTILS_H

#include "PrimesPolRho.h"
#include <Rcpp.h>
#include <gmp.h>

template <typename typeReturn>
void getPrimeFactors(std::int64_t &t, std::vector<typeReturn> &factors);

// 2^63 -->> 922337203685477580
constexpr double my63Max = 922337203685477580.0;

// std::sqrt(std::pow(2, 63)) -->> 3037000499.97605
constexpr std::int64_t Sqrt63Max = static_cast<std::int64_t>(3037000499.97);

/* Number of Miller-Rabin tests to run when not proving primality. */
constexpr int MR_REPS = 25;

int IsPrime(std::int64_t n);
void PollardRhoMpzT(mpz_t n, std::size_t a, std::vector<double> &factors);
void PollardRho(std::int64_t n, std::int64_t a, std::vector<int>& factors);

#endif
