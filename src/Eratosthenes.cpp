#include "NumbersUtils/PrimesSegSieve.h"
#include "NumbersUtils/Wheel.h"
#include <algorithm>
#include <vector>
#include <thread>
#include <cmath>
#include <deque>

// "PrimeSieve" implements a simple segmented version of the Sieve of
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here:
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve

namespace PrimeSieve {

    constexpr int L1_CACHE_SIZE = 32768;

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

    constexpr std::array<double, 15> PERCINC = {{0.2500, 0.1160, 0.1030, 0.0850, 0.0712,
                                                 0.0614, 0.0538, 0.0495, 0.0480, 0.0431,
                                                 0.0392, 0.0360, 0.0332, 0.0309, 0.0292}};

    constexpr std::array<double, 15> CUTPOINTS = {{40000.0,            120000.0,
                                                   1000000.0,          10000000.0,
                                                   100000000.0,        1000000000.0,
                                                   5000000000.0,       10000000000.0,
                                                   100000000000.0,     1000000000000.0,
                                                   10000000000000.0,   100000000000000.0,
                                                   1000000000000000.0, 5000000000000000.0,
                                                   10000000000000000.0}};

    // The following function is based off of the prime number theorem
    std::size_t EstimatePiPrime(double minNum, double maxNum) {
        const auto it = std::upper_bound(CUTPOINTS.cbegin(), CUTPOINTS.cend(), maxNum);
        const std::size_t myIndex = it - CUTPOINTS.cbegin();
        double dblRes = std::ceil((maxNum / std::log(maxNum)) * (1 + PERCINC[myIndex]));

        if (minNum > 1000) {
            dblRes -= std::floor((minNum / std::log(minNum)) * (1 + PERCINC[myIndex]));
        }

        return dblRes;
    }

    template <typename T>
    void PrimeSieveSmall(const std::vector<int> &sievePrimes,
                         std::vector<T> &primes, int minNum, int maxNum) {

        primes.reserve(EstimatePiPrime(
            static_cast<double>(minNum), static_cast<double>(maxNum)
        ));

        if (maxNum <= smallPrimeBase.back()) {
            int ind = 0;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}

            for (; smallPrimeBase[ind] <= maxNum; ++ind) {
                primes.push_back(smallPrimeBase[ind]);
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
                    primes.push_back(smallPrimeBase[ind]);
                }
            }

            int p = 1;
            int sqrPrime = 9;
            int myNum = 1 + lowerBnd;

            std::vector<int> nextStrt;
            nextStrt.reserve(sievePrimes.size());

            if (minNum > 2) {
                const int upperBnd = std::min(lowerBnd + segSize, maxNum);

                for (int myIndex = 0; sqrPrime <= upperBnd; ++p) {
                    if (lowerBnd > sqrPrime) {
                        const int remTest = lowerBnd % sievePrimes[p - 1];
                        if (remTest == 0) {
                            myIndex = sievePrimes[p - 1];
                        } else {
                            myIndex = sievePrimes[p - 1] - remTest;

                            if ((myIndex % 2) == 0) {
                                myIndex += sievePrimes[p - 1];
                            }
                        }
                    } else {
                        myIndex = sqrPrime - lowerBnd;
                    }

                    nextStrt.push_back(myIndex);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }

                for (std::size_t i = strt2310; i < nextStrt.size(); ++i) {
                    int j = nextStrt[i];

                    for (int k = sievePrimes[i] * 2; j < segSize; j += k) {
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
                            primes.push_back(myNum);
                        }
                    }

                    ++strtChunk;
                }

                if (upperBnd < flrMaxNum) {
                    int sieveInd = myNum - lowerBnd;

                    for (int q = strtChunk; q < N_WHEELS2310_PER_SEG; ++q) {
                        for (int w = 0; w < SZ_WHEEL2310; sieveInd += ARR_WHEEL2310[w + 7], w += 8) {
                            if (sieve[sieveInd]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 1]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 2]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 3]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 4]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 5]]) primes.push_back(sieveInd + lowerBnd);
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 6]]) primes.push_back(sieveInd + lowerBnd);
                        }
                    }

                    myNum = sieveInd + lowerBnd;
                } else {
                    for (int q = strtChunk; q < N_WHEELS2310_PER_SEG && myNum <= maxNum; ++q) {
                        for (int w = 0; w < SZ_WHEEL2310 && myNum <= maxNum;
                             myNum += ARR_WHEEL2310[w], ++w) {
                            if (sieve[myNum - lowerBnd]) {
                                primes.push_back(myNum);
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
                        for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                             j < segSize; j += k) {

                            sieve[j] = 0;
                        }

                        nextStrt[i] = (sievePrimes[i] * 2) - (
                            (segSize - nextStrt[i]) % (sievePrimes[i] * 2)
                        );
                    }
                } else {
                    for (int i = strt2310; i < halfSegSize2310Idx; ++i) {
                        for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                             j < segSize; j += k) {

                            sieve[j] = 0;
                        }

                        nextStrt[i] = (sievePrimes[i] * 2) - (
                            (segSize - nextStrt[i]) % (sievePrimes[i] * 2)
                        );
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
                    for (int w = 0; w < SZ_WHEEL2310; idx += ARR_WHEEL2310[w + 7], w += 8) {
                        if (sieve[idx]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 1]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 2]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 3]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 4]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 5]]) primes.push_back(idx + lowerBnd);
                        if (sieve[idx += ARR_WHEEL2310[w + 6]]) primes.push_back(idx + lowerBnd);
                    }
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

                for (int i = strt2310, size = nextStrt.size(); i < size; ++i) {
                    for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                         j < segSize; j += k) {

                        sieve[j] = 0;
                    }
                }

                for (int q = 0; q < N_WHEELS2310_PER_SEG && myNum <= maxNum; ++q) {
                    for (int w = 0; w < SZ_WHEEL2310 && myNum <= maxNum; myNum += ARR_WHEEL2310[w], ++w) {
                        if (sieve[myNum - lowerBnd]) {
                            primes.push_back(myNum);
                        }
                    }
                }
            }
        }
    }

    template <typename T>
    void PrimeSieveMedium(const std::vector<int> &sievePrimes,
                          std::vector<T> &primes, std::int_fast64_t minNum,
                          std::int_fast64_t maxNum) {

        constexpr std::size_t szWheel30030  = SZ_WHEEL30030;
        constexpr std::int_fast64_t sz30030 = NUM30030;
        const std::size_t nWheels = static_cast<std::size_t>(std::max(
            1.0, std::ceil(std::sqrt(static_cast<double>(maxNum)) / sz30030))
        );
        const std::int_fast64_t segSize = nWheels * sz30030;

        // We have to add primes.size() because we could have started in PrimeSieveSmall
        primes.reserve(
            EstimatePiPrime(static_cast<double>(minNum),
                            static_cast<double>(maxNum)) + primes.size()
        );

        const std::int_fast64_t flrMaxNum = segSize * (maxNum / segSize);
        std::vector<bool> sieve(segSize, true);

        std::int_fast64_t lowerBnd = segSize * (minNum / segSize);
        const std::int_fast64_t upperBnd = std::min(lowerBnd + segSize,
                                                    maxNum);
        std::int_fast64_t myNum = 1 + lowerBnd;

        std::size_t p = 1u;
        std::int_fast64_t sqrPrime = 9;
        std::vector<int> nextStrt;

        for (std::int_fast64_t myIndex = 0; sqrPrime <= upperBnd; ++p) {
            if (lowerBnd > sqrPrime) {
                const int remTest = lowerBnd % sievePrimes[p - 1u];

                if (remTest == 0) {
                    myIndex = sievePrimes[p - 1u];
                } else {
                    myIndex = sievePrimes[p - 1u] - remTest;
                    if ((myIndex % 2) == 0) {myIndex += sievePrimes[p - 1u];}
                }
            } else {
                myIndex = sqrPrime - lowerBnd;
            }

            nextStrt.push_back(myIndex);
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

        std::size_t strtChunk = 0u;

        if (minNum > myNum) {
            strtChunk = (minNum - myNum) / sz30030;
            myNum += strtChunk * sz30030;

            for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                 myNum += ARR_WHEEL30030[w], ++w) {

                if (myNum >= minNum && sieve[myNum - lowerBnd]) primes.push_back(myNum);
            }

            ++strtChunk;
        }

        if (upperBnd < flrMaxNum) {
            std::size_t sieveInd = myNum - lowerBnd;

            for (std::size_t q = strtChunk; q < nWheels; ++q) {
                for (auto w: ARR_WHEEL30030) {
                    if (sieve[sieveInd]) primes.push_back(sieveInd + lowerBnd);
                    sieveInd += w;
                }
            }

            myNum = sieveInd + lowerBnd;
        } else {
            for (std::size_t q = strtChunk; q < nWheels && myNum <= maxNum; ++q) {
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                    myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) primes.push_back(myNum);
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
                for (auto w: ARR_WHEEL30030) {
                    if (sieve[sieveInd]) primes.push_back(sieveInd + lowerBnd);
                    sieveInd += w;
                }
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

            for (std::size_t i = strt30030; i < nextStrt.size(); ++i) {
                for (int k = sievePrimes[i] * 2, j = nextStrt[i];
                     j < segSize; j += k) {
                    sieve[j] = false;
                }
            }

            for (std::size_t q = 0u; q < nWheels && myNum <= maxNum; ++q) {
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                     myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) primes.push_back(myNum);
                }
            }
        }
    }

    template <typename T>
    void PrimeSieveBig(const std::vector<int> &svPriOne,
                       const std::vector<int> &svPriTwo,
                       const std::vector<char> &check30030,
                       std::vector<T> &primes, std::size_t nBigSegs,
                       std::int_fast64_t minNum, std::int_fast64_t maxNum) {

        const int sz30030 = NUM30030;
        const int segSize = nBigSegs * sz30030;
        const std::size_t numWheelSegs = nBigSegs;
        const std::size_t szWheel30030 = SZ_WHEEL30030;

        // We have to add primes.size() because we could have started in PrimeSieveMedium
        primes.reserve(
            EstimatePiPrime(static_cast<double>(minNum),
                            static_cast<double>(maxNum)) + primes.size()
        );

        std::int_fast64_t lowerBnd = static_cast<std::int_fast64_t>(
            segSize * (minNum / segSize)
        );

        const std::size_t svPriOneSize = svPriOne.size();
        std::vector<int> nextStrt(svPriOneSize);

        for (std::size_t i = 0, idx = 0; i < svPriOneSize; ++i) {
            int remTest = lowerBnd % svPriOne[i];

            if (remTest == 0) {
                idx = svPriOne[i];
            } else {
                idx = svPriOne[i] - remTest;
                if (idx % 2 == 0) {idx += svPriOne[i];}
            }

            nextStrt[i] = idx;
        }

        const int myRange = (maxNum - lowerBnd) + 1;
        int numCacheSegs = myRange / segSize;

        double wholeTest = static_cast<double>(myRange) / segSize;
        if (wholeTest != static_cast<int>(wholeTest)) ++numCacheSegs;

        // Keeps track of which primes will be used in each interval
        std::deque<std::vector<int>> myBuckets(numCacheSegs,
                                               std::vector<int>());

        for (int i = 0, idx = 0, maxIndex = myRange + 1,
             size = svPriTwo.size(); i < size; ++i) {

            int remTest = (lowerBnd % svPriTwo[i]);

            if (remTest == 0) {
                idx = svPriTwo[i];
            } else {
                idx = (svPriTwo[i] - remTest);
                if ((idx % 2) == 0) {idx += svPriTwo[i];}
            }

            remTest = (idx % sz30030) - 1;
            int timesTwo = (2 * svPriTwo[i]);
            int remPrime = (timesTwo % sz30030);
            bool bKeepGoing = (idx <= maxIndex);

            // Populate rest of the buckets
            while (bKeepGoing) {
                while (check30030[remTest] && bKeepGoing) {
                    idx += timesTwo;
                    bKeepGoing = (idx <= maxIndex);
                    remTest += remPrime;
                    if (remTest >= sz30030) {remTest -= sz30030;}
                }

                int tempInd = idx;

                if (bKeepGoing) {
                    int divTest = (idx / segSize);
                    idx -= (divTest * segSize);
                    myBuckets[divTest].push_back(idx);
                }

                idx = tempInd + timesTwo;
                bKeepGoing = (idx <= maxIndex);
                remTest += remPrime;
                if (remTest >= sz30030) remTest -= sz30030;
            }
        }

        const std::int_fast64_t flrMaxNum = segSize * (maxNum / segSize);
        std::int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        std::int_fast64_t myNum = 1 + lowerBnd;

        // vector used for sieving
        std::vector<bool> sieve(segSize, true);
        std::size_t strt = 0u;

        if (minNum != lowerBnd) {
            for (std::size_t i = strt30030; i < svPriOneSize; ++i) {
                int j = nextStrt[i];

                for (int k = svPriOne[i] * 2; j < segSize; j += k) {
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

                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w) {
                    if (myNum >= minNum && sieve[myNum - lowerBnd]) primes.push_back(myNum);
                }

                ++strtChunk;
            }

            if (upperBnd < flrMaxNum) {
                std::int_fast64_t sieveInd = myNum - lowerBnd;

                for (std::size_t q = strtChunk; q < numWheelSegs; ++q) {
                    for (auto w: ARR_WHEEL30030) {
                        if (sieve[sieveInd]) primes.push_back(sieveInd + lowerBnd);
                        sieveInd += w;
                    }
                }

                myNum = sieveInd + lowerBnd;
            } else {
                for (std::size_t q = strtChunk; q < numWheelSegs && myNum <= maxNum; ++q) {
                    for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum;
                         myNum += ARR_WHEEL30030[w], ++w) {
                        if (sieve[myNum - lowerBnd]) primes.push_back(myNum);
                    }
                }
            }

            std::fill(sieve.begin(), sieve.end(), true);
            lowerBnd += segSize;
            ++strt;
        }

        for (std::size_t v = strt; lowerBnd < flrMaxNum; ++v, lowerBnd += segSize) {
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
                for (auto w: ARR_WHEEL30030) {
                    if (sieve[sieveInd]) primes.push_back(sieveInd + lowerBnd);
                    sieveInd += w;
                }
            }

            myNum = sieveInd + lowerBnd;
            std::fill(sieve.begin(), sieve.end(), true);
        }

        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (std::size_t i = strt30030; i < svPriOneSize; ++i) {
                for (int k = svPriOne[i] * 2, j = nextStrt[i]; j < segSize; j += k) {
                    sieve[j] = false;
                }
            }

            for (const auto &elem: myBuckets[0]) {
                sieve[elem] = false;
            }

            for (std::size_t q = 0u; (myNum <= maxNum) && (q < numWheelSegs); ++q) {
                for (std::size_t w = 0u; (myNum <= maxNum) && (w < szWheel30030); myNum += ARR_WHEEL30030[w], ++w) {
                    if (sieve[myNum - lowerBnd]) primes.push_back(myNum);
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
            if (bAddZero) sievePrimes.push_back(0);
            std::size_t ind = (bAddTwo) ? 0u : 1u;

            for (; smallPrimeBase[ind] <= sqrtBound; ++ind) {
                sievePrimes.push_back(smallPrimeBase[ind]);
            }

            if (bAddExtraPrime)
                sievePrimes.push_back(smallPrimeBase[ind]);
        } else {
            int sqrtSqrtBound = static_cast<int>(std::sqrt(static_cast<double>(sqrtBound)));
            std::vector<int> sqrtSievePrimes;
            sqrtSmallPrimes(sqrtSqrtBound, sqrtSievePrimes);

            // The number, 225, comes from the observation that the largest prime
            // gap less than 100 million is 219 @ 47,326,693. This is important
            // because we need to guarantee that we obtain the smallest prime
            // greater than sqrt(n). We know that the largest number that this
            // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266

            int myLower = 3;
            int myUpper = bAddExtraPrime ? sqrtBound + 225 : sqrtBound;
            const std::size_t sqrtReserve = EstimatePiPrime(1.0, static_cast<double>(myUpper));
            sievePrimes.reserve(sqrtReserve);

            if (bAddZero) sievePrimes.push_back(0);
            if (bAddTwo) myLower = 1;
            PrimeSieveSmall(sqrtSievePrimes, sievePrimes, myLower, myUpper);
        }
    }

    template <typename T>
    void PrimeWorker(const std::vector<int> &svPriMain,
                     const std::vector<int> &svPriOne,
                     const std::vector<int> &svPriTwo,
                     const std::vector<char> &check30030,
                     std::vector<T> &primes, std::int_fast64_t minNum,
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

    template <typename T>
    void PrimeSieveMain(std::vector<std::vector<T>> &primeList,
                        std::vector<T> &primes, std::int_fast64_t minNum,
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
        std::int_fast64_t sqrtMedCut = maxCores * (1024 * 1024) * 6;

        const int sqrtBound = std::sqrt(static_cast<double>(maxNum));
        std::size_t nBigSegs = 1u;

        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) {nThreads = maxThreads;}
            sqrtMedCut /= nThreads;
        }

        std::vector<int> svPriOne;
        std::vector<int> svPriTwo;
        std::vector<int> svPriMain;

        const std::int_fast64_t medCut = sqrtMedCut * sqrtMedCut;
        std::vector<char> check30030(NUM30030, 1);
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
            nBigSegs = ((maxThreads / nThreads) * 2 * 1024 * 1024) /  L1_CACHE_SIZE;

            std::size_t ind = 0u;
            const int limitOne = static_cast<int>(nBigSegs * segUnitSize);
            const std::size_t svMainSize = svPriMain.size();

            // Get the primes that are guaranteed to mark an index in every segment interval
            for (; ind < svMainSize && (2 * svPriMain[ind]) < limitOne; ++ind) {
                svPriOne.push_back(svPriMain[ind]);
            }

            // Get the rest
            for (; ind < svMainSize; ++ind) {
                svPriTwo.push_back(svPriMain[ind]);
            }

            for (std::size_t i = 0u, ind = 0u; i < SZ_WHEEL30030;
                 ind += ARR_WHEEL30030[i], ++i) {
                check30030[ind] = 0;
            }

        } else if (maxNum > smallCut) {
            nBigSegs = static_cast<std::size_t>(std::ceil(
                static_cast<double>(sqrtBound) / segUnitSize));
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
        if (maxThreads < 2) Parallel = false;

        if (Parallel) {
            if (testNumThreads > 1) {
                if (testNumThreads < nThreads) {nThreads = testNumThreads;}
            } else {
                Parallel = false;
            }
        }

        if (Parallel) {
            std::vector<std::thread> threads;
            std::int_fast64_t upperBnd;
            std::int_fast64_t lowerBnd = minNum;

            const std::int_fast64_t chunkSize = segSize * std::floor(
                static_cast<double>(myRange / nThreads) / segSize
            );
            upperBnd = segSize * std::floor(minNum / segSize) + chunkSize;

            for (int idx = 0; idx < (nThreads - 1); lowerBnd = upperBnd, upperBnd += chunkSize, ++idx) {
                threads.emplace_back(std::cref(PrimeWorker<T>), std::cref(svPriMain),
                                     std::cref(svPriOne), std::cref(svPriTwo),
                                     std::cref(check30030), std::ref(primeList[idx]),
                                     lowerBnd, upperBnd, smallCut, medCut, nBigSegs);
            }

            threads.emplace_back(std::cref(PrimeWorker<T>), std::cref(svPriMain),
                                 std::cref(svPriOne), std::cref(svPriTwo),
                                 std::cref(check30030), std::ref(primeList.back()),
                                 lowerBnd, maxNum, smallCut, medCut, nBigSegs);

            for (auto& thr: threads) {
                thr.join();
            }
        } else {
            PrimeWorker(svPriMain, svPriOne, svPriTwo, check30030,
                        primes, minNum, maxNum, smallCut, medCut, nBigSegs);
        }
    }
}

template void PrimeSieve::PrimeSieveMain(std::vector<std::vector<int>>&,
                                         std::vector<int>&, std::int_fast64_t,
                                         std::int_fast64_t, bool&, int, int, int);
template void PrimeSieve::PrimeSieveMain(std::vector<std::vector<double>>&,
                                         std::vector<double>&, std::int_fast64_t,
                                         std::int_fast64_t, bool&, int, int, int);
template void PrimeSieve::PrimeSieveMain(std::vector<std::vector<std::int64_t>>&,
                                         std::vector<std::int64_t>&, std::int_fast64_t,
                                         std::int_fast64_t, bool&, int, int, int);

template void PrimeSieve::sqrtBigPrimes(int, bool, bool, bool, std::vector<int>&);
template void PrimeSieve::sqrtBigPrimes(int, bool, bool, bool, std::vector<double>&);
template void PrimeSieve::sqrtBigPrimes(int, bool, bool, bool, std::vector<std::int64_t>&);
