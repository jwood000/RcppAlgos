#include "NumbersUtils/PrimeSieveCount.h"
#include "NumbersUtils/PrimesSegSieve.h"
#include "NumbersUtils/PhiTinyLookup.h"
#include "NumbersUtils/Eratosthenes.h"
#include "CleanConvert.h"
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>
#include <cmath>

// "MainPrimeCount" based off of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

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

    constexpr int MAX_A = 100;
    std::array<std::vector<std::uint16_t>, MAX_A> phiCache;
    std::mutex theBlocker;

    // Increment MAX_A, so we can have easier access to indexing. E.g.
    // when a = 3 (i.e. the number of primes is 3), instead of accessing
    // the third entry in phiTiny like phiTiny[a - 1], we simply use "a"
    // as our index (i.e. phiTiny[a]). The first entry in phiTiny will
    // be empty. The same goes for phiPrimes and phiPi.
    std::vector<std::int64_t> phiPrimes;
    std::vector<std::int64_t> phiPi;

    // The arrays below are utilized by phiTiny
    // primeProds[n] = \prod_{i=1}^{n} primes[i]
    // myTotients[n] = \prod_{i=1}^{n} (primes[i] - 1)
    std::int64_t phiTinyCalc(std::int64_t x, std::int64_t a) {
        static constexpr std::array<int, 7> myTotients = {{1, 1, 2, 8, 48, 480, 5760}};
        static constexpr std::array<int, 7> primeProds = {{1, 2, 6, 30, 210, 2310, 30030}};

        std::int64_t pp = primeProds[a];
        return (x / pp) * myTotients[a] + phiTiny[a][x % pp];
    }

    void updateCache(std::uint64_t x, std::uint64_t a, std::int64_t mySum) {
        if (a < phiCache.size() &&
            x <= std::numeric_limits<std::uint16_t>::max()) {
            // Protect phiCache while its being updated
            std::lock_guard<std::mutex> guard(theBlocker);
            if (x >= phiCache[a].size()) phiCache[a].resize(x + 1, 0);
            phiCache[a][x] = static_cast<std::uint16_t>(std::abs(mySum));
        }
    }

    bool isCached(std::uint64_t x, std::uint64_t a) {
        return a < phiCache.size() &&
               x < phiCache[a].size() &&
               phiCache[a][x];
    }

    bool isPix(std::int64_t x, std::int64_t a) {
        return x < static_cast<std::int64_t>(phiPi.size()) &&
               x < (phiPrimes[a + 1] * phiPrimes[a + 1]);
    }

    std::int64_t getStrt(std::int64_t y) {
        static constexpr std::array<int, 7> myTinyPrimes = {{0, 2, 3, 5, 7, 11, 13}};
        static constexpr std::array<int, 13> myTinyPi = {{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5}};

        if (y >= myTinyPrimes.back())
            return phiTinySize;
        else
            return myTinyPi[y];
    }

    template <std::int64_t SIGN>
    std::int64_t phiWorker(std::int64_t x, std::int64_t a) {
        if (x <= phiPrimes[a]) {
            return SIGN;
        } else if (a <= phiTinySize) {
            return phiTinyCalc(x, a) * SIGN;
        } else if (isPix(x, a)) {
            return (phiPi[x] - a + 1) * SIGN;
        } else if (isCached(x, a)) {
            return phiCache[a][x] * SIGN;
        } else {
            std::size_t sqrtx = static_cast<std::size_t>(std::sqrt(static_cast<double>(x)));
            std::int64_t piSqrtx = a;
            std::int64_t strt = getStrt(sqrtx);

            if (sqrtx < phiPi.size()) {
                piSqrtx = std::min(static_cast<std::int64_t>(phiPi[sqrtx]), a);
            }

            std::int64_t mySum = (piSqrtx - a) * SIGN + phiTinyCalc(x, strt) * SIGN;

            for (std::int64_t i = strt; i < piSqrtx; ++i) {
                std::int64_t x2 = x / phiPrimes[i + 1];

                if (isPix(x2, i)) {
                    mySum += (phiPi[x2] - i + 1) * -SIGN;
                } else {
                    mySum += phiWorker<-SIGN>(x2, i);
                }
            }

            updateCache(x, a, mySum);
            return mySum;
        }
    }

    void phiForeman(std::int64_t &mySum, std::int64_t lowerBound,
                    std::int64_t upperBound, std::int64_t x) {

        for (std::int64_t i = lowerBound; i < upperBound; ++i) {
            mySum += phiWorker<-1>(x / phiPrimes[i + 1], i);
        }
    }

    const double getChunkFactor(std::int64_t x) {
        static constexpr std::array<double, 9> nums = {{1e10, 1e12, 2e13, 5e13, 8e13, 1e14, 5e14, 1e15, 1e16}};
        static constexpr std::array<double, 9> factor = {{1.3, 1.2, 1.1, 1.07, 1.05, 1.01, 1.007, 1.006, 1.005}};

        const auto it = std::upper_bound(nums.cbegin(), nums.cend(), static_cast<double>(x));
        return std::log(factor[it - nums.cbegin()]);
    }

    std::int64_t phiMain(std::int64_t x, std::int64_t a, int nThreads, bool Parallel) {

        const std::int64_t sqrtx = static_cast<std::int64_t>(std::sqrt(static_cast<double>(x)));
        const std::int64_t piSqrtx = std::min(static_cast<std::int64_t>(phiPi[sqrtx]), a);
        const std::int64_t strt = getStrt(sqrtx);
        std::int64_t mySum = phiTinyCalc(x, strt) + piSqrtx - a;

        if (Parallel) {
            constexpr double divLim = 10000000000.0;    //  1e10

            // We know x, t, and we need to solve for n.
            // x / (1.5^t)^n < 1e10 -->>> 1e10 * (1.5^t)^n = x --->>> x / 1e10 = (1.5^t)^n
            // log(x / 1e10) = n * log(1.5^t) --->> log(x / 1e10) / log(1.5^t) :
            //                 log(x / 1e10) / (t * log(1.5))
            std::int64_t myRange = (piSqrtx - strt) + 1;
            std::int64_t lower = strt;
            std::int64_t upper;

            if (x > divLim) {
                const double chunkFactor = getChunkFactor(x);
                const int nLoops = 1 + std::ceil(std::log(x / divLim) / (nThreads * chunkFactor));
                const double multOne = std::exp(std::log(myRange) / (nThreads * nLoops));
                const int firstStep = std::ceil(std::log10(Significand53 / x)) + 1;
                std::int64_t base = lower;

                //*********************** Begin of firstThreads *****************************
                // Typically, the size of multOne will be very small resulting in not very
                // meaningful chunk sizes in the first few loops. Because of this, we
                // gurantee that each of the first threads will have at least firstStep
                // iterations and that the threads following will have at least firstStep + 1

                // Here we have m = multOne, p = power, and s = firstStep:
                // m^(p +  1) - m^p = m^p * (m - 1) > (s + 1)  -->>  m^p > (s + 1) / (m - 1)  -->>
                // p * log(m) > log((s + 1) / (m - 1))  -->>  p > log((s + 1) / (m - 1)) / log(m)
                double power = std::ceil(std::log((firstStep + 1) / (multOne - 1)) / std::log(multOne)) + 1;

                // Same as above: lower + (upper - lower) / s   and given:   m^p + lower = upper
                //   -->>  (lower + m^p / s)  -->>  firstLoops = ceil(((m^p + lower) / s) / nThreads)
                const int firstLoops = std::ceil(((std::pow(multOne, power) +
                                                 base) / firstStep) / nThreads);
                upper = lower + firstStep;
                std::vector<std::int64_t> mySums;

                for (int i = 0; i < firstLoops; ++i) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    std::vector<std::thread> threads;

                    for (int j = 0; j < nThreads; lower = upper, upper += firstStep, ++j) {
                        threads.emplace_back(phiForeman, std::ref(intermediate[j]),
                                             lower, upper, x);
                    }

                    for (auto& thr: threads) {
                        thr.join();
                    }

                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                //*********************** End of firstThreads *****************************
                //*************************************************************************
                //*********************** Begin of midThreads *****************************

                while (std::pow(multOne, power) < (upper - base))
                    ++power;

                double dblChunk = std::pow(multOne, power);
                upper = static_cast<std::int64_t>(dblChunk) + base;

                while ((dblChunk * std::pow(multOne, nThreads - 1) + base) < piSqrtx) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    std::vector<std::thread> threads;

                    for (int j = 0; j < nThreads; lower = upper, ++j,
                            dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                        threads.emplace_back(phiForeman, std::ref(intermediate[j]),
                                             lower, upper, x);
                    }

                    for (auto& thr: threads) {
                        thr.join();
                    }

                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                std::vector<std::int64_t> intermediate(nThreads, 0);
                std::vector<std::thread> threads;

                for (int j = 0; j < (nThreads - 1) && upper < piSqrtx; lower = upper, ++j,
                        dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                    threads.emplace_back(phiForeman, std::ref(intermediate[j]),
                                         lower, upper, x);
                }

                threads.emplace_back(phiForeman, std::ref(intermediate.back()),
                                     lower, piSqrtx, x);

                for (auto& thr: threads) {
                    thr.join();
                }

                mySums.insert(mySums.end(), intermediate.cbegin(),
                              intermediate.cend());

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            } else {
                std::vector<std::int64_t> mySums(nThreads, 0);
                std::vector<std::thread> threads;
                std::int64_t chunk = myRange / nThreads;
                upper = lower + chunk - 1;

                for (int j = 0; j < (nThreads - 1); lower = upper, upper += chunk, ++j) {
                    threads.emplace_back(phiForeman, std::ref(mySums[j]),
                                         lower, upper, x);
                }

                threads.emplace_back(phiForeman, std::ref(mySums.back()),
                                     lower, piSqrtx, x);

                for (auto& thr: threads) {
                    thr.join();
                }

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            }
        } else {
            phiForeman(mySum, strt, piSqrtx, x);
        }

        return mySum;
    }

    // All values verified by Kim Walisch's primecount library (nThreads = 8)
    //  10^9 -->>          50,847,534   -->>   711.5 microseconds
    // 10^10 -->>         455,052,511   -->>   3.25  milliseconds
    // 10^11 -->>       4,118,054,813   -->>   12.45 milliseconds
    // 10^12 -->>      37,607,912,018   -->>   67.95 milliseconds
    // 10^13 -->>     346,065,536,839   -->>   507.4 milliseconds
    // 10^14 -->>   3,204,941,750,802   -->>   3.867 seconds
    // 10^15 -->>  29,844,570,422,669   -->>  28.407 seconds
    // MAX VALUE (2^53 - 1) -->>
    //            252,252,704,148,404   -->> 213.575 seconds

    std::int64_t MainPrimeCount(std::int64_t n, int nThreads = 1, int maxThreads = 1) {

        const std::int64_t sqrtBound = std::sqrt(static_cast<double>(n));
        std::vector<std::int64_t> resetPhiPrimes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, true, false, true, resetPhiPrimes);
        phiPrimes = resetPhiPrimes;

        phiPi.resize(sqrtBound + 1);
        std::int64_t count = 0;
        const std::int64_t maxPrime = phiPrimes.back();

        for (std::int64_t i = 1; i <= maxPrime; ++i) {
            if (i >= phiPrimes[count + 1]) {
                ++count;
            }

            phiPi[i] = count;
        }

        for (std::int64_t i = (maxPrime + 1); i <= sqrtBound; ++i)
            phiPi[i] = count;

        bool Parallel = false;

        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) nThreads = maxThreads;
            if ((maxThreads < 2) || (n < 1e7)) Parallel = false;
        }

        const std::int64_t piSqrt = PiPrime(sqrtBound);
        const std::int64_t phiSqrt = phiMain(n, piSqrt, nThreads, Parallel);
        const std::int64_t int64result = piSqrt + phiSqrt - 1;

        return int64result;
    }
}

SEXP PrimeCountCpp(SEXP Rn, SEXP RNumThreads, SEXP RmaxThreads) {
    double dblNum;
    CleanConvert::convertPrimitive(Rn, dblNum, VecType::Numeric, "n");
    const std::int64_t n = static_cast<std::int64_t>(dblNum);

    int nThreads = 1;
    int maxThreads = 1;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    if (n < 100000) {
        if (n < 10) {
            if (n == 1)
                return Rf_ScalarInteger(0);
            else if (n == 2)
                return Rf_ScalarInteger(1);
            else if (n < 5)
                return Rf_ScalarInteger(2);
            else if (n < 7)
                return Rf_ScalarInteger(3);
            else
                return Rf_ScalarInteger(4);
        }

        return Rf_ScalarInteger(static_cast<int>(PrimeCounting::PiPrime(n)));
    }

    if (!Rf_isNull(RNumThreads)) {
        CleanConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    std::int64_t result = PrimeCounting::MainPrimeCount(n, nThreads, maxThreads);

    if (result > std::numeric_limits<int>::max()) {
        return Rf_ScalarReal(static_cast<double>(result));
    } else {
        return Rf_ScalarInteger(static_cast<int>(result));
    }
}
