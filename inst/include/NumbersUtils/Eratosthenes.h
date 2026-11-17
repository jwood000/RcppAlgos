#pragma once

#include <cmath>
#include <vector>
#include <cstdint>

namespace PrimeSieve{
    std::size_t EstimatePiPrime(double minNum, double maxNum);

    template <typename T>
    void PrimeSieveMain(std::vector<std::vector<T>> &primeList,
                        std::vector<T> &primes, std::int_fast64_t minNum,
                        std::int_fast64_t maxNum, bool &Parallel,
                        int nThreads = 1, int maxThreads = 1,
                        int maxCores = 1);

    // Internal inclusive interval counter using the same Small/Medium/Big
    // sieve and worker policy as generation. Caller validates 1 <= minNum <=
    // maxNum <= 2^53 - 1 and bounds maxNum - minNum to at most 1e9.
    // Caller supplies 1 <= nThreads <= maxThreads and maxCores >= 1.
    // Stores one count per useful worker, never generated prime values.
    std::int64_t CountPrimes(std::int_fast64_t minNum,
                             std::int_fast64_t maxNum, int nThreads = 1,
                             int maxThreads = 1, int maxCores = 1);

    template <typename T>
    void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                       bool bAddTwo, std::vector<T> &sievePrimes);
}
