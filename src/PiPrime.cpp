#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Eratosthenes.h"
#include "NumbersUtils/PrimesSegSieve.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <stdexcept>

namespace PrimeCounting {
    namespace {

        // Big prepares bucket entries for each window. Bound that working storage
        // on long ranges without retaining a vector of the primes found.
        constexpr std::int64_t SIEVE_WINDOW = 1000000000;

    }

    std::int64_t PiPrimeRange(std::int64_t minNum, std::int64_t maxNum,
                              int nThreads, int maxThreads, int maxCores) {

        if (minNum < 0 || maxNum < 0 ||
            minNum > Significand53Int || maxNum > Significand53Int) {
            throw std::out_of_range(
                    "PiPrime bounds must be between 0 and 2^53 - 1"
            );
        }

        if (minNum > maxNum || maxNum < 2) {
            return 0;
        }

        if (maxNum <= smallPrimeBase.back()) {
            const auto first = std::lower_bound(
                smallPrimeBase.cbegin(), smallPrimeBase.cend(), minNum
            );

            const auto last = std::upper_bound(
                first, smallPrimeBase.cend(), maxNum
            );

            return last - first;
        }

        minNum = std::max<std::int64_t>(2, minNum);
        maxThreads = std::max(1, maxThreads);
        nThreads = std::max(1, std::min(nThreads, maxThreads));
        maxCores = std::max(1, maxCores);

        std::int64_t count = 0;

        while (minNum <= maxNum) {
            const std::int64_t upper = std::min(
                maxNum, minNum + SIEVE_WINDOW - 1
            );

            count += PrimeSieve::CountPrimes(
                minNum, upper, nThreads, maxThreads, maxCores
            );

            minNum = upper + 1;
        }

        return count;
    }

    std::int64_t PiPrime(std::int64_t maxNum) {
        return PiPrimeRange(0, maxNum);
    }
}
