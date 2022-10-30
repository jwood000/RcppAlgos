#pragma once

#include <vector>
#include <cstdint>

namespace MotleyPrimes {

    template <typename T, typename U>
    void EulerPhiSieve(T m, U retN, T offsetStrt,
                       const std::vector<T> &primes,
                       std::vector<T> &numSeq,
                       U* EulerPhis);
}
