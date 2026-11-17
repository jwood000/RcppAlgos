#include "NumbersUtils/PrimesSegSieve.h"
#include "NumbersUtils/Wheel.h"
#include "NumbersUtils/ThreadBatch.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

// "PrimeSieve" implements a simple segmented version of the Sieve of
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here:
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve

namespace PrimeSieve {

    //************************ Prime Counting Estimates ************************
    // These numbers were obtained empirically using the prime number theorem
    // along with a prime counting function. Values were computed in each
    // range below in cutPoints. For example in the range [40000, 120000),
    // if we increment by 10 and calculate the following:
    //                       pi(x) - x / log(x) = e
    // We then calculate the percentage of x / log(x) or : e / (x / log(x))
    // for every value and find the max:
    //             a <- sapply(seq(40000, 120000, 10), primeCount)
    //             b <- sapply(seq(40000, 120000, 10), function(x) x/log(x))
    //             myDiff <- a - b;  max(myDiff / b) ## [1] 0.1153694

    constexpr std::array<double, 15> PERCINC = {
        {0.2500, 0.1160, 0.1030, 0.0850, 0.0712,
         0.0614, 0.0538, 0.0495, 0.0480, 0.0431,
         0.0392, 0.0360, 0.0332, 0.0309, 0.0292}
    };

    constexpr std::array<double, 15> CUTPOINTS = {
        {40000.0,            120000.0,
         1000000.0,          10000000.0,
         100000000.0,        1000000000.0,
         5000000000.0,       10000000000.0,
         100000000000.0,     1000000000000.0,
         10000000000000.0,   100000000000000.0,
         1000000000000000.0, 5000000000000000.0,
         10000000000000000.0}
    };

    // The following function is based off of the prime number theorem
    std::size_t EstimatePiPrime(double minNum, double maxNum) {
        const auto it = std::upper_bound(
            CUTPOINTS.cbegin(), CUTPOINTS.cend(), maxNum
        );

        const std::size_t idx = it - CUTPOINTS.cbegin();
        double dblRes = std::ceil((maxNum / std::log(maxNum)) * (1 + PERCINC[idx]));

        if (minNum > 1000) {
            dblRes -= std::floor((minNum / std::log(minNum)) * (1 + PERCINC[idx]));
        }

        return dblRes;
    }

    // Output overloads keep marking identical for generation and counting.
    // Counting has no reserve estimate, compaction buffer, or prime vector.
    template <typename T>
    void ReservePrimes(std::vector<T>& primes, std::int_fast64_t minNum,
                       std::int_fast64_t maxNum, bool append = false) {
        primes.reserve(EstimatePiPrime(static_cast<double>(minNum),
                                      static_cast<double>(maxNum)) +
                       (append ? primes.size() : 0));
    }

    inline void ReservePrimes(std::int64_t&, std::int_fast64_t,
                              std::int_fast64_t, bool = false) {}

    template <typename T>
    void AppendPrime(std::vector<T>& primes, std::int_fast64_t value) {
        primes.push_back(value);
    }

    inline void AppendPrime(std::int64_t& count, std::int_fast64_t) {
        ++count;
    }

    template <typename T, typename Sieve, typename Index, std::size_t N>
    void AppendWheel(std::vector<T>& primes, const Sieve& sieve, Index& index,
                      std::int_fast64_t lowerBnd, const char (&wheel)[N]) {
        // Rejected candidates are overwritten; only survivors are appended.
        // One slot per candidate bounds all writes, including rejected ones.
        T compact[N];
        T* writeTo = compact;

        for (auto gap : wheel) {
            *writeTo = index + lowerBnd;
            writeTo += (sieve[index] != 0);
            index += gap;
        }

        primes.insert(primes.end(), compact, writeTo);
    }

    template <typename Sieve, typename Index, std::size_t N>
    void AppendWheel(std::int64_t& count, const Sieve& sieve, Index& index,
                      std::int_fast64_t, const char (&wheel)[N]) {

        for (auto gap : wheel) {
            count += (sieve[index] != 0);
            index += gap;
        }
    }

    template <typename Output>
    void PrimeSieveSmall(const std::vector<int> &sievePrimes,
                         Output &primes, int minNum, int maxNum) {

        ReservePrimes(primes, minNum, maxNum);

        if (maxNum <= smallPrimeBase.back()) {
            int ind = 0;
            const int smPriBaseSize = smallPrimeBase.size();

            for (; ind < smPriBaseSize && smallPrimeBase[ind] < minNum; ++ind) {}

            for (; ind < smPriBaseSize && smallPrimeBase[ind] <= maxNum; ++ind) {
                AppendPrime(primes, smallPrimeBase[ind]);
            }
        } else {
            constexpr int segSize = Almost2310L1Cache;
            const int flrMaxNum   = segSize * (maxNum / segSize);

            int lowerBnd = segSize * (minNum / segSize);
            std::vector<char> sieve(segSize, 1);

            if (minNum < wheel2310PrimeLim) {
                sieve[1] = 0;
                std::size_t ind = 0u;
                for (; smallPrimeBase[ind] < minNum; ++ind) {}

                for (; smallPrimeBase[ind] < wheel2310PrimeLim; ++ind) {
                    AppendPrime(primes, smallPrimeBase[ind]);
                }
            }

            int p = 1;
            int sqrPrime = 9;
            int myNum = 1 + lowerBnd;

            std::vector<int> nextStrt;
            nextStrt.reserve(sievePrimes.size());

            if (minNum > 2) {
                const int upperBnd = std::min(lowerBnd + segSize, maxNum);

                for (int idx = 0; sqrPrime <= upperBnd; ++p) {
                    if (minNum > sqrPrime) {
                        const int prime = sievePrimes[p - 1u];
                        const int remTest = minNum % prime;
                        int first = minNum + (remTest ? prime - remTest : 0);

                        if (first % 2 == 0) {
                            first += prime;
                        }

                        idx = first - lowerBnd;
                    } else {
                        idx = sqrPrime - lowerBnd;
                    }

                    nextStrt.push_back(idx);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }

                const std::int_fast64_t offsetLim =
                    std::min(segSize, maxNum - lowerBnd + 1);

                for (std::size_t i = strt2310; i < nextStrt.size(); ++i) {
                    int j = nextStrt[i];

                    for (int k = sievePrimes[i] * 2; j < offsetLim; j += k) {
                        sieve[j] = 0;
                    }

                    nextStrt[i] = j - segSize;
                }

                int strtChunk = 0;

                if (minNum > myNum) {
                    strtChunk = (minNum - myNum) / NUM2310;
                    myNum += strtChunk * NUM2310;

                    for (int w = 0; w < SZ_WHEEL2310 && myNum <= maxNum;
                         myNum += ARR_WHEEL2310[w], ++w) {
                        if (myNum >= minNum && sieve[myNum - lowerBnd]) {
                            AppendPrime(primes, myNum);
                        }
                    }

                    ++strtChunk;
                }

                if (upperBnd < flrMaxNum) {
                    int sieveInd = myNum - lowerBnd;

                    for (int q = strtChunk; q < N_WHEELS2310_PER_SEG; ++q) {

                        AppendWheel(primes, sieve, sieveInd, lowerBnd, ARR_WHEEL2310);
                    }

                    myNum = sieveInd + lowerBnd;
                } else {
                    for (int q = strtChunk; q < N_WHEELS2310_PER_SEG && myNum <= maxNum; ++q) {
                        for (int w = 0; w < SZ_WHEEL2310 && myNum <= maxNum;
                             myNum += ARR_WHEEL2310[w], ++w) {
                            if (sieve[myNum - lowerBnd]) {
                                AppendPrime(primes, myNum);
                            }
                        }
                    }
                }

                std::fill(sieve.begin(), sieve.end(), 1);
                lowerBnd += segSize;
            }

            for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
                for (; sqrPrime <= (lowerBnd + segSize); ++p) {
                    nextStrt.push_back(sqrPrime - lowerBnd);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }

                if (lowerBnd + segSize < (segSize * segSize / 4)) {
                    for (std::size_t i = strt2310; i < nextStrt.size(); ++i) {
                        int j = nextStrt[i];

                        for (int k = sievePrimes[i] * 2; j < segSize; j += k) {
                            sieve[j] = 0;
                        }

                        nextStrt[i] = j - segSize;
                    }
                } else {
                    for (int i = strt2310; i < halfSegSize2310Idx; ++i) {
                        int j = nextStrt[i];

                        for (int k = sievePrimes[i] * 2; j < segSize; j += k) {
                            sieve[j] = 0;
                        }

                        nextStrt[i] = j - segSize;
                    }

                    for (std::size_t i = halfSegSize2310Idx; i < nextStrt.size(); ++i) {
                        if (nextStrt[i] > segSize) {
                            nextStrt[i] -= segSize;
                        } else {
                            sieve[nextStrt[i]] = 0;
                            nextStrt[i] = (nextStrt[i] + sievePrimes[i] * 2) - segSize;
                        }
                    }
                }

                for (int q = 0, idx = myNum - lowerBnd; q < N_WHEELS2310_PER_SEG; ++q) {
                    AppendWheel(primes, sieve, idx, lowerBnd, ARR_WHEEL2310);
                }

                myNum += segSize;
                std::fill(sieve.begin(), sieve.end(), 1);
            }

            // Get remaining primes that are greater than flrMaxNum and less than maxNum
            if (lowerBnd < maxNum) {
                for (; sqrPrime <= maxNum; ++p) {
                    nextStrt.push_back(sqrPrime - lowerBnd);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }

                const std::int_fast64_t offsetLim = maxNum - lowerBnd + 1;

                for (int i = strt2310, size = nextStrt.size(); i < size; ++i) {
                    for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                         j < offsetLim; j += k) {

                        sieve[j] = 0;
                    }
                }

                for (int q = 0; q < N_WHEELS2310_PER_SEG && myNum <= maxNum; ++q) {
                    for (int w = 0; w < SZ_WHEEL2310 && myNum <= maxNum;
                         myNum += ARR_WHEEL2310[w], ++w) {
                        if (sieve[myNum - lowerBnd]) {
                            AppendPrime(primes, myNum);
                        }
                    }
                }
            }
        }
    }

    template <typename Output>
    void PrimeSieveMedium(const std::vector<int> &sievePrimes,
                          Output &primes, std::int_fast64_t minNum,
                          std::int_fast64_t maxNum) {

        constexpr std::size_t szWheel30030  = SZ_WHEEL30030;
        constexpr std::int_fast64_t sz30030 = NUM30030;
        const std::size_t nWheels = static_cast<std::size_t>(std::max(
            1.0, std::ceil(std::sqrt(static_cast<double>(maxNum)) / sz30030))
        );
        const std::int_fast64_t segSize = nWheels * sz30030;

        // Earlier Small/Medium phases may already have emitted primes.
        ReservePrimes(primes, minNum, maxNum, true);

        const std::int_fast64_t flrMaxNum = segSize * (maxNum / segSize);
        std::vector<bool> sieve(segSize, true);

        std::int_fast64_t lowerBnd = segSize * (minNum / segSize);
        const std::int_fast64_t upperBnd = std::min(
            lowerBnd + segSize, maxNum
        );
        std::int_fast64_t myNum = 1 + lowerBnd;

        std::size_t p = 1u;
        std::int_fast64_t sqrPrime = 9;
        std::vector<int> nextStrt;

        for (std::int_fast64_t idx = 0; sqrPrime <= upperBnd; ++p) {
            if (minNum > sqrPrime) {
                const int prime = sievePrimes[p - 1u];
                const int remTest = minNum % prime;

                std::int_fast64_t first =
                    minNum + (remTest ? prime - remTest : 0);

                if (first % 2 == 0) {
                    first += prime;
                }

                idx = first - lowerBnd;
            } else {
                idx = sqrPrime - lowerBnd;
            }

            nextStrt.push_back(idx);
            sqrPrime = static_cast<std::int_fast64_t>(sievePrimes[p]) *
                static_cast<std::int_fast64_t>(sievePrimes[p]);
        }

        const std::int_fast64_t offsetLim =
            std::min(segSize, maxNum - lowerBnd + 1);

        for (std::size_t i = strt30030; i < nextStrt.size(); ++i) {
            int j = nextStrt[i];

            for (int k = sievePrimes[i] * 2; j < offsetLim; j += k) {
                sieve[j] = false;
            }

            nextStrt[i] = j - segSize;
        }

        std::size_t strtChunk = 0u;

        if (minNum > myNum) {
            strtChunk = (minNum - myNum) / sz30030;
            myNum += strtChunk * sz30030;

            for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                 myNum += ARR_WHEEL30030[w], ++w) {
                if (myNum >= minNum && sieve[myNum - lowerBnd]) {
                    AppendPrime(primes, myNum);
                }
            }

            ++strtChunk;
        }

        if (upperBnd < flrMaxNum) {
            std::size_t sieveInd = myNum - lowerBnd;

            for (std::size_t q = strtChunk; q < nWheels; ++q) {

                AppendWheel(primes, sieve, sieveInd, lowerBnd, ARR_WHEEL30030);
            }

            myNum = sieveInd + lowerBnd;
        } else {
            for (std::size_t q = strtChunk; q < nWheels && myNum <= maxNum; ++q) {
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                    myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) AppendPrime(primes, myNum);
                }
            }
        }

        std::fill(sieve.begin(), sieve.end(), true);
        lowerBnd += segSize;

        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            for (; sqrPrime <= (lowerBnd + segSize); ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = static_cast<std::int_fast64_t>(sievePrimes[p]) *
                    static_cast<std::int_fast64_t>(sievePrimes[p]);
            }

            for (std::size_t i = strt30030; i < nextStrt.size(); ++i) {
                int j = nextStrt[i];

                for (int k = sievePrimes[i] * 2; j < segSize; j += k) {
                    sieve[j] = false;
                }

                nextStrt[i] = j - segSize;
            }

            for (std::size_t q = 0u, sieveInd = myNum - lowerBnd; q < nWheels; ++q) {
                AppendWheel(primes, sieve, sieveInd, lowerBnd, ARR_WHEEL30030);
            }

            myNum += segSize;
            std::fill(sieve.begin(), sieve.end(), true);
        }

        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = static_cast<std::int_fast64_t>(sievePrimes[p]) *
                    static_cast<std::int_fast64_t>(sievePrimes[p]);
            }

            const std::int_fast64_t offsetLim = maxNum - lowerBnd + 1;

            for (std::size_t i = strt30030; i < nextStrt.size(); ++i) {
                for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                     j < offsetLim; j += k) {
                    sieve[j] = false;
                }
            }

            for (std::size_t q = 0u; q < nWheels && myNum <= maxNum; ++q) {
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                     myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) AppendPrime(primes, myNum);
                }
            }
        }
    }

    template <typename Output>
    void PrimeSieveBig(const std::vector<int> &svPriOne,
                       const std::vector<int> &svPriTwo,
                       const std::vector<char> &check30030,
                       Output &primes, std::size_t nBigSegs,
                       std::int_fast64_t minNum, std::int_fast64_t maxNum) {

        const int sz30030 = NUM30030;
        const int segSize = static_cast<int>(nBigSegs * sz30030);
        const std::size_t numWheelSegs = nBigSegs;
        const std::size_t szWheel30030 = SZ_WHEEL30030;

        // Earlier Small/Medium phases may already have emitted primes.
        ReservePrimes(primes, minNum, maxNum, true);

        std::int_fast64_t lowerBnd = static_cast<std::int_fast64_t>(
            segSize * (minNum / segSize)
        );

        const std::size_t svPriOneSize = svPriOne.size();
        std::vector<int> nextStrt(svPriOneSize);

        for (std::size_t i = 0, idx = 0; i < svPriOneSize; ++i) {
            int remTest = minNum % svPriOne[i];
            std::int_fast64_t first =
                minNum + (remTest ? svPriOne[i] - remTest : 0);

            if (first % 2 == 0) {
                first += svPriOne[i];
            }

            idx = first - lowerBnd;
            nextStrt[i] = idx;
        }

        // The aligned prefix can make this wider than the requested interval.
        // Only positions within one segment are stored as int offsets.
        const std::int_fast64_t myRange = (maxNum - lowerBnd) + 1;
        const std::int_fast64_t cacheSegments =
            myRange / segSize + (myRange % segSize != 0);
        const std::size_t numCacheSegs = static_cast<std::size_t>(cacheSegments);
        if (static_cast<std::int_fast64_t>(numCacheSegs) != cacheSegments) {
            throw std::length_error("Prime sieve bucket count exceeds size_t");
        }

        // Keeps track of which primes will be used in each interval
        std::deque<std::vector<int>> myBuckets(numCacheSegs,
                                               std::vector<int>());

        for (std::size_t i = 0; i < svPriTwo.size(); ++i) {

            int remTest = minNum % svPriTwo[i];
            std::int_fast64_t first =
                minNum + (remTest ? svPriTwo[i] - remTest : 0);

            if (first % 2 == 0) {
                first += svPriTwo[i];
            }

            std::int_fast64_t idx = first - lowerBnd;
            remTest = (idx % sz30030) - 1;
            int timesTwo = (2 * svPriTwo[i]);
            int remPrime = (timesTwo % sz30030);
            bool bKeepGoing = (idx < myRange);

            // Populate rest of the buckets
            while (bKeepGoing) {
                while (check30030[remTest] && bKeepGoing) {
                    idx += timesTwo;
                    bKeepGoing = (idx < myRange);
                    remTest += remPrime;

                    if (remTest >= sz30030) {
                        remTest -= sz30030;
                    }
                }

                const std::int_fast64_t tempInd = idx;

                if (bKeepGoing) {
                    const std::size_t divTest = static_cast<std::size_t>(idx / segSize);
                    myBuckets[divTest].push_back(static_cast<int>(idx % segSize));
                }

                idx = tempInd + timesTwo;
                bKeepGoing = (idx < myRange);
                remTest += remPrime;

                if (remTest >= sz30030) {
                    remTest -= sz30030;
                }
            }
        }

        const std::int_fast64_t flrMaxNum = segSize * (maxNum / segSize);
        std::int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        std::int_fast64_t myNum = 1 + lowerBnd;

        // vector used for sieving
        std::vector<bool> sieve(segSize, true);

        if (minNum != lowerBnd) {
            const int offsetLim = static_cast<int>(
                std::min<std::int_fast64_t>(
                    segSize, maxNum - lowerBnd + 1
                )
            );

            for (std::size_t i = strt30030; i < svPriOneSize; ++i) {
                int j = nextStrt[i];

                for (int k = svPriOne[i] * 2; j < offsetLim; j += k) {
                    sieve[j] = false;
                }

                nextStrt[i] = j - segSize;
            }

            for (const auto &elem: myBuckets[0]) {
                sieve[elem] = false;
            }

            myBuckets.pop_front();
            std::size_t strtChunk = 0u;

            if (minNum > myNum) {
                strtChunk = (minNum - myNum) / NUM30030;
                myNum += strtChunk * NUM30030;

                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                     myNum += ARR_WHEEL30030[w], ++w) {
                    if (myNum >= minNum && sieve[myNum - lowerBnd]) {
                        AppendPrime(primes, myNum);
                    }
                }

                ++strtChunk;
            }

            if (upperBnd < flrMaxNum) {
                std::int_fast64_t sieveInd = myNum - lowerBnd;

                for (std::size_t q = strtChunk; q < numWheelSegs; ++q) {
                    AppendWheel(primes, sieve, sieveInd, lowerBnd, ARR_WHEEL30030);
                }

                myNum = sieveInd + lowerBnd;
            } else {
                for (std::size_t q = strtChunk; q < numWheelSegs && myNum <= maxNum; ++q) {
                    for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                         myNum += ARR_WHEEL30030[w], ++w) {
                        if (sieve[myNum - lowerBnd]) AppendPrime(primes, myNum);
                    }
                }
            }

            std::fill(sieve.begin(), sieve.end(), true);
            lowerBnd += segSize;
        }

        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            for (std::size_t i = strt30030; i < svPriOneSize; ++i) {
                int j = nextStrt[i];

                for (int k = (svPriOne[i] * 2); j < segSize; j += k) {
                    sieve[j] = false;
                }

                nextStrt[i] = j - segSize;
            }

            for (const auto &elem: myBuckets[0]) {
                sieve[elem] = false;
            }

            myBuckets.pop_front();
            std::int_fast64_t sieveInd = myNum - lowerBnd;

            for (std::size_t q = 0u; q < numWheelSegs; ++q) {

                AppendWheel(primes, sieve, sieveInd, lowerBnd, ARR_WHEEL30030);
            }

            myNum = sieveInd + lowerBnd;
            std::fill(sieve.begin(), sieve.end(), true);
        }

        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            const int offsetLim = static_cast<int>(maxNum - lowerBnd + 1);

            for (std::size_t i = strt30030; i < svPriOneSize; ++i) {
                for (int k = svPriOne[i] * 2, j = nextStrt[i];
                     j < offsetLim; j += k) {
                    sieve[j] = false;
                }
            }

            for (const auto &elem: myBuckets[0]) {
                sieve[elem] = false;
            }

            for (std::size_t q = 0u; (myNum <= maxNum) && (q < numWheelSegs); ++q) {
                for (std::size_t w = 0u; (myNum <= maxNum) && (w < szWheel30030); myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) AppendPrime(primes, myNum);
                }
            }
        }
    }

    inline void sqrtSmallPrimes(int sqrtBound, std::vector<int> &sievePrimes) {
        int ind = 1;

        for (; smallPrimeBase[ind] <= sqrtBound; ++ind) {
            sievePrimes.push_back(smallPrimeBase[ind]);
        }

        sievePrimes.push_back(smallPrimeBase[ind]);
    }

    template <typename T>
    void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                       bool bAddTwo, std::vector<T> &sievePrimes) {

        if (sqrtBound < smallPrimeBase.back()) {
            const auto first = smallPrimeBase.cbegin() + (bAddTwo ? 0 : 1);
            auto last = std::upper_bound(
                first, smallPrimeBase.cend(), sqrtBound
            );

            if (bAddExtraPrime) {
                ++last;
            }

            sievePrimes.reserve(
                sievePrimes.size() + (last - first) + (bAddZero ? 1u : 0u)
            );

            if (bAddZero) {
                sievePrimes.push_back(0);
            }

            sievePrimes.insert(sievePrimes.end(), first, last);
        } else {
            // The number, 225, comes from the fact that the largest prime
            // gap less than 100 million is 219 @ 47,326,693. This is important
            // because we need to guarantee that we obtain the smallest prime
            // greater than sqrt(n). We know that the largest number that this
            // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266

            int myUpper = bAddExtraPrime ? sqrtBound + 225 : sqrtBound;
            int sqrtSqrtBound =
                static_cast<int>(std::sqrt(static_cast<double>(myUpper)));
            std::vector<int> sqrtSievePrimes;
            sqrtSmallPrimes(sqrtSqrtBound, sqrtSievePrimes);

            int myLower = 3;
            const std::size_t sqrtReserve = EstimatePiPrime(
                1.0, static_cast<double>(myUpper)
            );
            sievePrimes.reserve(sqrtReserve);

            if (bAddZero) sievePrimes.push_back(0);
            if (bAddTwo) myLower = 1;
            PrimeSieveSmall(sqrtSievePrimes, sievePrimes, myLower, myUpper);
        }
    }

    template <typename Output>
    void PrimeWorker(const std::vector<int> &svPriMain,
                     const std::vector<int> &svPriOne,
                     const std::vector<int> &svPriTwo,
                     const std::vector<char> &check30030,
                     Output &primes, std::int_fast64_t minNum,
                     std::int_fast64_t maxNum, std::int_fast64_t smallCut,
                     std::int_fast64_t medCut, std::size_t nBigSegs) {

        if (maxNum > medCut) {
            // PrimeSieveBig is not meant for small values, so we
            // first get them with PrimeSieveMed/Small
            if (minNum < smallCut) {
                PrimeSieveSmall(svPriMain, primes, minNum, smallCut);
                PrimeSieveMedium(svPriMain, primes, smallCut, medCut);
                PrimeSieveBig(svPriOne, svPriTwo, check30030, primes, nBigSegs, medCut, maxNum);
            } else if (minNum < medCut) {
                PrimeSieveMedium(svPriMain, primes, minNum, medCut);
                PrimeSieveBig(svPriOne, svPriTwo, check30030, primes, nBigSegs, medCut, maxNum);
            } else {
                PrimeSieveBig(svPriOne, svPriTwo, check30030, primes, nBigSegs, minNum, maxNum);
            }
        } else if (maxNum > smallCut) {
            if (minNum < smallCut) {
                PrimeSieveSmall(svPriMain, primes, minNum, smallCut);
                PrimeSieveMedium(svPriMain, primes, smallCut, maxNum);
            } else {
                PrimeSieveMedium(svPriMain, primes, minNum, maxNum);
            }
        } else {
            PrimeSieveSmall(svPriMain, primes, minNum, maxNum);
        }
    }

    template <typename Output>
    void SieveMain(std::vector<Output> &primeList,
                        Output &primes, std::int_fast64_t minNum,
                        std::int_fast64_t maxNum, bool &Parallel,
                        int nThreads, int maxThreads, int maxCores) {

        const std::int_fast64_t myRange = maxNum - minNum;
        std::int_fast64_t smallCut = Almost2310L1Cache - 100;
        smallCut *= smallCut;

        // Here, we estimate the maximum sieving interval. Based off of empirical
        // evidence we start seeing improved performance in PrimeSieveBig after
        // (sqrtMedCut)^2 vs PrimeSieveMedium. We suspect this has to do with the
        // amount of L3 Cache available. On our main testing device, we have 6MiB.
        // If you were to calculate the sieving size using the formula below, we
        // would obtain (4 * 1024^2 * 6) / 8 (std::vector<bool> possible memory
        // optimization of packing elements into individual bits) ~= 3MiB which
        // gives us room to store other important objects for sieving.
        std::int_fast64_t sqrtMedCut =
            static_cast<std::int_fast64_t>(maxCores) * (1024 * 1024) * 6;

        const int sqrtBound = std::sqrt(static_cast<double>(maxNum));
        std::size_t nBigSegs = 1u;

        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;

            if (nThreads > maxThreads) {
                nThreads = maxThreads;
            }

            sqrtMedCut /= nThreads;
        }

        std::vector<int> svPriOne;
        std::vector<int> svPriTwo;
        std::vector<int> svPriMain;

        // Above sqrt(2^53 - 1), every supported input uses Small/Medium.
        // Avoid squaring a hardware-derived value that can exceed int64.
        // PrimeWorker requires ordered cutoffs, even when a native caller's
        // thread budget is much larger than its physical-core budget.
        const std::int_fast64_t medCut = std::max<std::int_fast64_t>(
            smallCut,
            sqrtMedCut > Sqrt53Max ? Significand53Int : sqrtMedCut * sqrtMedCut
        );

        std::vector<char> check30030(maxNum > medCut ? NUM30030 : 0, 1);
        std::int_fast64_t segUnitSize = NUM30030;

        // bool bAddZero, bool bAddExtraPrime, bool bAddTwo
        sqrtBigPrimes(sqrtBound, false, true, false, svPriMain);

        if (maxNum > medCut) {
            // Similar to sqrtMedCut above (see comment). This formulation
            // was determined empirically. The idea is that as the number
            // of threads increases, we want to ensure that we don't pollute
            // the cache memory. Each thread will need enough space for all
            // of the data structures necessary to produce primes, thus
            // we reduce the segment size as nThreads increases.
            const std::int_fast64_t requestedSegs =
                (static_cast<std::int_fast64_t>(maxThreads / nThreads) *
                 2 * 1024 * 1024) / L1_CACHE_SIZE;

            // Dense marking advances once past the segment. Keep that final
            // position, and each initial nextStrt, representable as int.
            const std::int_fast64_t maxSegment =
                std::numeric_limits<int>::max() - 2LL * svPriMain.back();

            nBigSegs = std::min(
                requestedSegs,
                maxSegment / segUnitSize
            );

            const std::int_fast64_t limitOne = nBigSegs * segUnitSize;

            // limitOne is even: 2 * prime < limitOne means prime < limitOne / 2.
            const auto split = std::lower_bound(
                svPriMain.cbegin(), svPriMain.cend(), limitOne / 2
            );

            svPriOne.assign(svPriMain.cbegin(), split);
            svPriTwo.assign(split, svPriMain.cend());

            for (std::size_t i = 0u, idx = 0u; i < SZ_WHEEL30030;
                 idx += ARR_WHEEL30030[i], ++i) {
                check30030[idx] = 0;
            }

        } else if (maxNum > smallCut) {
            nBigSegs = static_cast<std::size_t>(
                std::ceil(static_cast<double>(sqrtBound) / segUnitSize)
            );
        } else {
            segUnitSize = NUM2310;
            nBigSegs = N_WHEELS2310_PER_SEG;
        }

        nBigSegs = (nBigSegs < 1) ? 1u : nBigSegs;
        const std::int_fast64_t segSize = nBigSegs * segUnitSize;

        // There is some overhead for setting up the data structures
        // for each thread, so we ensure each thread has enough work
        // to outweigh adding another thread. Determined empirically.
        // The tests that were performed showed that one thread was
        // just as efficient as two threads when myRange was equal
        // to 4 * segSize. This is especially true for PrimeSieveBig.
        const std::int_fast64_t testNumThreads = myRange / (segSize * 4);

        if (maxThreads < 2) {
            Parallel = false;
        }

        if (Parallel && testNumThreads > 1 && testNumThreads < nThreads) {
            nThreads = testNumThreads;
        } else if (Parallel && testNumThreads <= 1) {
            Parallel = false;
        }

        if (Parallel) {
            // Count-only callers allocate after useful-worker reduction.
            // Generation retains its caller-owned slots and final back slot.
            if (primeList.size() < static_cast<std::size_t>(nThreads)) {
                primeList.resize(nThreads);
            }

            // Shared sieve inputs and caller-owned output vectors outlive the
            // batch. It joins every launched worker before propagating either
            // a worker exception or a later thread-creation failure.
            PrimeCounting::ThreadBatch batch(nThreads);
            std::int_fast64_t upperBnd;
            std::int_fast64_t lowerBnd = minNum;

            const std::int_fast64_t chunkSize = segSize * std::floor(
                static_cast<double>(myRange / nThreads) / segSize
            );

            upperBnd = segSize * std::floor(minNum / segSize) + chunkSize;

            for (int idx = 0; idx < (nThreads - 1);
                 lowerBnd = upperBnd, upperBnd += chunkSize, ++idx) {

                batch.launch([&, idx, lowerBnd, upperBnd]() {
                    PrimeWorker(svPriMain, svPriOne, svPriTwo, check30030,
                                primeList[idx], lowerBnd, upperBnd,
                                smallCut, medCut, nBigSegs);
                });
            }

            batch.launch([&, lowerBnd]() {
                PrimeWorker(svPriMain, svPriOne, svPriTwo, check30030,
                            primeList.back(), lowerBnd, maxNum,
                            smallCut, medCut, nBigSegs);
            });

            batch.finish();
        } else {
            PrimeWorker(svPriMain, svPriOne, svPriTwo, check30030,
                        primes, minNum, maxNum, smallCut, medCut, nBigSegs);
        }
    }

    template <typename T>
    void PrimeSieveMain(std::vector<std::vector<T>> &primeList,
                        std::vector<T> &primes, std::int_fast64_t minNum,
                        std::int_fast64_t maxNum, bool &Parallel,
                        int nThreads, int maxThreads, int maxCores) {

        SieveMain(primeList, primes, minNum, maxNum, Parallel,
                  nThreads, maxThreads, maxCores);
    }

    // Caller supplies 1 <= nThreads <= maxThreads and maxCores >= 1.
    std::int64_t CountPrimes(std::int_fast64_t minNum,
                             std::int_fast64_t maxNum, int nThreads,
                             int maxThreads, int maxCores) {

        std::vector<std::int64_t> counts;
        std::int64_t serialCount = 0;
        bool parallel = false;

        SieveMain(counts, serialCount, minNum, maxNum, parallel,
                  nThreads, maxThreads, maxCores);

        if (!parallel) {
            return serialCount;
        }

        std::int64_t result = 0;

        for (auto count : counts) {
            result += count;
        }

        return result;
    }
}

template void PrimeSieve::PrimeSieveMain(
    std::vector<std::vector<int>>&, std::vector<int>&, std::int_fast64_t,
    std::int_fast64_t, bool&, int, int, int
);
template void PrimeSieve::PrimeSieveMain(
    std::vector<std::vector<double>>&, std::vector<double>&, std::int_fast64_t,
    std::int_fast64_t, bool&, int, int, int
);
template void PrimeSieve::PrimeSieveMain(
    std::vector<std::vector<std::int64_t>>&, std::vector<std::int64_t>&,
    std::int_fast64_t, std::int_fast64_t, bool&, int, int, int
);

template void PrimeSieve::sqrtBigPrimes(
    int, bool, bool, bool, std::vector<int>&
);
template void PrimeSieve::sqrtBigPrimes(
    int, bool, bool, bool, std::vector<double>&
);
template void PrimeSieve::sqrtBigPrimes(
    int, bool, bool, bool, std::vector<std::int64_t>&
);
