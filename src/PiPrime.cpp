#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/PrimesSegSieve.h"
#include "NumbersUtils/Wheel.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace PrimeCounting {

    // PiPrime is very similar to the PrimeSieveSmall only we are not
    // considering a range. That is, we are only concerned with finding
    // the number of primes less than maxNum. We are also only counting
    // primes instead of generating them.
    std::int64_t PiPrime(std::int64_t maxNum) {

        constexpr int segSize = Almost210L1Cache;
        constexpr int nWheels = N_WHEELS210_PER_SEG;
        constexpr int szWheel210 = SZ_WHEEL210;
        const int sqrtBound = std::sqrt(static_cast<double>(maxNum));

        // the wheel already has the first 4 primes marked as
        // false, so we need to account for them here. N.B.
        // the calling functions has checks for cases where
        // maxNum is less than 11, so no need to check here.
        std::int64_t count = 4;

        std::vector<int> smallPrimes;
        std::vector<int> nextStrt;
        const int flrMaxNum = segSize * (maxNum / segSize);

        std::size_t ind = 1;

        for (; smallPrimeBase[ind] <= sqrtBound; ++ind) {
            smallPrimes.push_back(smallPrimeBase[ind]);
        }

        smallPrimes.push_back(smallPrimeBase[ind]);
        std::vector<char> sieve(segSize, 1);
        sieve[1] = 0;

        int sqrPrime = 9;
        int lowerBnd = 0;
        int myNum = 1;
        std::size_t p = 1;

        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            for (; sqrPrime <= (lowerBnd + segSize); ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }

            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                for (int k = smallPrimes[i] * 2, j = nextStrt[i];
                     j < segSize; j += k) {

                    sieve[j] = 0;
                }

                nextStrt[i] = (smallPrimes[i] * 2) - (
                    (segSize - nextStrt[i]) % (smallPrimes[i] * 2)
                );
            }

            for (int q = 0, idx = myNum - lowerBnd; q < nWheels; ++q) {
                for (auto w: ARR_WHEEL210) {
                    if (sieve[idx]) ++count;
                    idx += w;
                }
            }

            myNum += segSize;
            std::fill(sieve.begin(), sieve.end(), 1);
        }

        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }

            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                for (int k = smallPrimes[i] * 2, j = nextStrt[i];
                     j < segSize; j += k) {
                    sieve[j] = 0;
                }
            }

            for (int q = 0; q < nWheels && myNum <= maxNum; ++q) {
                for (int w = 0; w < szWheel210 && myNum <= maxNum;
                     myNum += ARR_WHEEL210[w], ++w) {
                    if (sieve[myNum - lowerBnd]) ++count;
                }
            }
        }

        return count;
    }

}
