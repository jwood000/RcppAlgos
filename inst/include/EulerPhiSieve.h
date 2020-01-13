#ifndef EULER_PHI_SIEVE_H
#define EULER_PHI_SIEVE_H

#include <vector>
#include <cstdint>

namespace MotleyPrimes {
    template <typename typeInt, typename typeReturn, typename typeRcpp>
    void EulerPhiSieve(typeInt m, typeReturn retN, typeInt offsetStrt,
                       const std::vector<typeInt> &primes,
                       std::vector<typeInt> &numSeq,
                       typeRcpp &EulerPhis);
}

#endif
