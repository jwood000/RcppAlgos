#ifndef PRIME_FACTORIZE_SIEVE_H
#define PRIME_FACTORIZE_SIEVE_H

#include <cstddef>
#include <vector>

namespace MotleyPrimes {
    template <typename typeInt>
    void PrimeFactorizationSieve(typeInt m, typeInt retN, typeInt offsetStrt,
                                 const std::vector<typeInt> &primes,
                                 std::vector<std::vector<typeInt>> &primeList);
}

#endif
