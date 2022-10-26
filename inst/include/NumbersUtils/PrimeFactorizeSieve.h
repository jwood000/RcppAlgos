#pragma once

#include <cstddef>
#include <vector>

namespace MotleyPrimes {
    template <typename T>
    void PrimeFactorizationSieve(T m, T retN, T offsetStrt,
                                 const std::vector<T> &primes,
                                 std::vector<std::vector<T>> &primeList);
}
