#ifndef ERATOSTHENES_H
#define ERATOSTHENES_H

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

    template <typename T>
    void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                       bool bAddTwo, std::vector<T> &sievePrimes);
}

#endif
