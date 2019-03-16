#include "MotleyPrimes.h"
#include "CleanConvert.h"
#include "PhiTinyLookup.h"
#include <RcppThread.h>
#include <mutex>
#include <array>

// "MasterPrimeCount" based off of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

namespace PrimeCounting {
    
    // PiPrime is very similar to the PrimeSieveSmall only we are not
    // considering a range. That is, we are only concerned with finding
    // the number of primes less than maxNum. We are also only counting
    // primes instead of generating them.
    int64_t PiPrime (int64_t maxNum) {
        
        const int_fast64_t segSize = Almost210L1Cache;
        const std::size_t nWheels = N_WHEELS210_PER_SEG;
        const std::size_t szWheel210 = SZ_WHEEL210;
        const int sqrtBound = static_cast<int>(std::sqrt(maxNum));
        
        // the wheel already has the first 4 primes marked as
        // false, so we need to account for them here. N.B.
        // the calling functions has checks for cases where
        // maxNum is less than 11, so no need to check here.
        int64_t count = 4;
        
        std::vector<int64_t> smallPrimes, nextStrt;
        const int64_t flrMaxNum = segSize * std::floor(maxNum / segSize);
        
        std::size_t ind = 1;
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            smallPrimes.push_back(smallPrimeBase[ind]);
        
        smallPrimes.push_back(smallPrimeBase[ind]);
        std::vector<char> sieve(segSize, 1);
        sieve[1] = 0;
        
        int64_t sqrPrime = 9;
        int64_t lowerBnd = 0, upperBnd = segSize, myNum = 1;
        unsigned long int p = 1;
        
        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            upperBnd = lowerBnd + segSize;
            
            for (; sqrPrime <= upperBnd; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int64_t j = nextStrt[i];
                for (int64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < nWheels; ++q)
                for (std::size_t w = 0; w < szWheel210; myNum += ARR_WHEEL210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        ++count;
                    
            std::fill(sieve.begin(), sieve.end(), 1);
        }
        
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int64_t j = nextStrt[i];
                for (int64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < nWheels && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += ARR_WHEEL210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        ++count;
        }
        
        return count;
    }

    constexpr int MAX_A = 100;
    std::array<std::vector<uint16_t>, MAX_A> phiCache;
    std::mutex theBlocker;
    
    // Increment MAX_A, so we can have easier access
    // to indexing.  E.g when a = 3 (i.e. the number)
    // of primes is 3), instead of accessing the third
    // entry in phiTiny like phiTiny[a - 1], we simply
    // use "a" as our index (i.e. phiTiny[a]). The
    // first entry in phiTiny will be empty. The same
    // goes for phiPrimes and phiPi.
    std::vector<int64_t> phiPrimes;
    std::vector<int64_t> phiPi;
    
    // The arrays below are utilized by phiTiny
    // primeProds[n] = \prod_{i=1}^{n} primes[i]
    // myTotients[n] = \prod_{i=1}^{n} (primes[i] - 1)
    constexpr std::array<int, 7> myTinyPrimes = {{0, 2, 3, 5, 7, 11, 13}};
    constexpr std::array<int, 7> primeProds = {{1, 2, 6, 30, 210, 2310, 30030}};
    constexpr std::array<int, 7> myTotients = {{1, 1, 2, 8, 48, 480, 5760}};
    constexpr std::array<int, 13> myTinyPi = {{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5}};
    
    int64_t phiTinyCalc(int64_t x, int64_t a) {
        int64_t pp = primeProds[a];
        return (x / pp) * myTotients[a] + phiTiny[a][x % pp];
    }
    
    void updateCache(uint64_t x, uint64_t a, int64_t mySum) {
        if (a < phiCache.size() &&
            x <= std::numeric_limits<uint16_t>::max()) {
            // Protect phiCache while its being updated
            std::lock_guard<std::mutex> guard(theBlocker);
            if (x >= phiCache[a].size()) phiCache[a].resize(x + 1, 0);
            phiCache[a][x] = static_cast<uint16_t>(std::abs(mySum));
        }
    }
    
    bool isCached(uint64_t x, uint64_t a) {
        return a < phiCache.size() &&
               x < phiCache[a].size() &&
               phiCache[a][x];
    }
    
    bool isPix(int64_t x, int64_t a) {
        return x < static_cast<int64_t>(phiPi.size()) &&
               x < (phiPrimes[a + 1] * phiPrimes[a + 1]);
    }
    
    int64_t getStrt(int64_t y) {
        if (y >= myTinyPrimes.back())
            return phiTinySize;
        else
            return myTinyPi[y];
    }
    
    template <int SIGN>
    int64_t phiWorker(int64_t x, int64_t a) {
        if (x <= phiPrimes[a]) {
            return SIGN;
        } else if (a <= phiTinySize) {
            return phiTinyCalc(x, a) * SIGN;
        } else if (isPix(x, a)) {
            return (phiPi[x] - a + 1) * SIGN;
        } else if (isCached(x, a)) {
            return phiCache[a][x] * SIGN;
        } else {
            int64_t sqrtx = static_cast<int64_t>(std::sqrt(x));
            int64_t piSqrtx = a;
            int64_t strt = getStrt(sqrtx);
            
            if (sqrtx < static_cast<int64_t>(phiPi.size()))
                piSqrtx = std::min(static_cast<int64_t>(phiPi[sqrtx]), a);
            
            int64_t mySum = 0;
            mySum += (piSqrtx - a) * SIGN;
            mySum += phiTinyCalc(x, strt) * SIGN;
            
            for (int64_t i = strt; i < piSqrtx; ++i) {
                int64_t x2 = x / phiPrimes[i + 1];
                
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
    
    int64_t phiSlave(int64_t lowerBound, int64_t upperBound, int64_t x) {
        
        int64_t mySum = 0;
        
        for (int64_t i = lowerBound; i < upperBound; ++i)
            mySum += phiWorker<-1>(x / phiPrimes[i + 1], i);
        
        return mySum;
    }
    
    const double getChunkFactor(int64_t x) {
        const std::vector<double> nums = {1e10, 1e12, 2e13, 5e13, 8e13, 1e14, 5e14, 1e15, 1e16};
        const std::vector<double> factor = {1.3, 1.2, 1.1, 1.07, 1.05, 1.01, 1.007, 1.006, 1.005};
        std::vector<double>::const_iterator it = std::upper_bound(nums.cbegin(),
                                                                  nums.cend(),
                                                                  static_cast<double>(x));
        return std::log(factor[it - nums.cbegin()]);
    }
    
    int64_t phiMaster(int64_t x, int64_t a, int nThreads, bool Parallel) {
        
        const int64_t sqrtx = static_cast<int64_t>(std::sqrt(x));
        const int64_t piSqrtx = std::min(static_cast<int64_t>(phiPi[sqrtx]), a);
        const int64_t strt = getStrt(sqrtx);
        int64_t mySum = phiTinyCalc(x, strt) + piSqrtx - a;
        
        if (Parallel) {
            const double divLim = 10000000000.0;    //  1e10
            
            // We know x, t, and we need to solve for n.
            // x / (1.5^t)^n < 1e10 -->>> 1e10 * (1.5^t)^n = x --->>> x / 1e10 = (1.5^t)^n
            // log(x / 1e10) = n * log(1.5^t) --->> log(x / 1e10) / log(1.5^t) :
            //                 log(x / 1e10) / (t * log(1.5))
            int64_t myRange = (piSqrtx - strt) + 1;
            int64_t lower = strt;
            int64_t upper, chunk;
            RcppThread::ThreadPool pool(nThreads);
            std::vector<std::future<int64_t>> myFutures;
            
            if (x > divLim) {
                const double chunkFactor = getChunkFactor(x);
                const int nLoops = 1 + std::ceil(std::log(x / divLim) / (nThreads * chunkFactor));
                const double multOne = std::exp(std::log(myRange) / (nThreads * nLoops));
                const int firstStep = std::ceil(std::log10(Significand53 / x)) + 1;
                int64_t base = lower;
                
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
                
                for (int i = 0; i < firstLoops; ++i) {
                    for (int j = 0; j < nThreads; lower = upper, upper += firstStep, ++j)
                        myFutures.push_back(pool.pushReturn(phiSlave, lower, upper, x));
                    
                    pool.wait();
                }
                
                //*********************** End of firstThreads *****************************
                //*************************************************************************
                //*********************** Begin of midThreads *****************************
                
                while (std::pow(multOne, power) < (upper - base))
                    ++power;
                
                double dblChunk = std::pow(multOne, power);
                upper = static_cast<int64_t>(dblChunk) + base;

                while ((dblChunk * std::pow(multOne, nThreads - 1) + base) < piSqrtx) {
                    for (int j = 0; j < nThreads; lower = upper, ++j, 
                            dblChunk *= multOne, upper = static_cast<int64_t>(dblChunk) + base) {
                        myFutures.push_back(pool.pushReturn(phiSlave, lower, upper, x));
                    }

                    pool.wait();
                }

                for (int j = 0; j < (nThreads - 1) && upper < piSqrtx; lower = upper, ++j, 
                        dblChunk *= multOne, upper = static_cast<int64_t>(dblChunk) + base) {
                    myFutures.push_back(pool.pushReturn(phiSlave, lower, upper, x));
                }

                myFutures.push_back(pool.pushReturn(phiSlave, lower, piSqrtx, x));

                for (std::size_t j = 0; j < myFutures.size(); ++j)
                    mySum += myFutures[j].get();

                pool.join();
                
            } else {
                chunk = myRange / nThreads;
                upper = lower + chunk - 1;
                
                for (int j = 0; j < (nThreads - 1); lower = upper, upper += chunk, ++j)
                    myFutures.push_back(pool.pushReturn(phiSlave, lower, upper, x));
                
                myFutures.push_back(pool.pushReturn(phiSlave, lower, piSqrtx, x));
                
                for (std::size_t j = 0; j < myFutures.size(); ++j)
                    mySum += myFutures[j].get();
                
                pool.join();
            }
        } else {
            mySum += phiSlave(strt, piSqrtx, x);
        }
        
        return mySum;
    }
    
    // All values verified by Kim Walisch's primecount library (nThreads = 8)
    //  10^9 -->>          50,847,534   -->>   811.6 microseconds
    // 10^10 -->>         455,052,511   -->>   5.362 milliseconds
    // 10^11 -->>       4,118,054,813   -->>   21.45 milliseconds
    // 10^12 -->>      37,607,912,018   -->>   119.4 milliseconds
    // 10^13 -->>     346,065,536,839   -->>   868.4 milliseconds
    // 10^14 -->>   3,204,941,750,802   -->>   6.724 seconds
    // 10^15 -->>  29,844,570,422,669   -->>  49.554 seconds
    // MAX VALUE (2^53 - 1) -->> 
    //            252,252,704,148,404   -->> 352.862 seconds
    
    int64_t MasterPrimeCount(int64_t n, int nThreads = 1, int maxThreads = 1) {
        
        const int64_t sqrtBound = static_cast<int64_t>(std::sqrt(n));
        std::vector<int64_t> resetPhiPrimes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, true, false, true, resetPhiPrimes);
        phiPrimes = resetPhiPrimes;
        
        phiPi.resize(sqrtBound + 1);
        int64_t count = 0;
        const int64_t maxPrime = phiPrimes.back();
        
        for (int64_t i = 1; i <= maxPrime; ++i) {
            if (i >= phiPrimes[count + 1])
                ++count;
            phiPi[i] = count;
        }
        
        for (int64_t i = (maxPrime + 1); i <= sqrtBound; ++i)
            phiPi[i] = count;
        
        bool Parallel = false;
        
        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) {nThreads = maxThreads;}
            if ((maxThreads < 2) || (n < 1e7)) {Parallel = false;}
        }
        
        const int64_t piSqrt = PiPrime(sqrtBound);
        const int64_t phiSqrt = phiMaster(n, piSqrt, nThreads, Parallel);
        const int64_t int64result = piSqrt + phiSqrt - 1;
        
        return int64result;
    }
}

//[[Rcpp::export]]
SEXP PrimeCountRcpp(SEXP Rn, SEXP RNumThreads, int maxThreads) {
    double dblNum;
    CleanConvert::convertPrimitive(Rn, dblNum, "n");

    if (dblNum < 1 || dblNum > Significand53)
        Rcpp::stop("n must be a positive number less than 2^53");

    const int64_t n = dblNum;
    
    if (n < 100000) {
        if (n < 10) {
            if (n == 1)
                return Rcpp::wrap(0);
            else if (n == 2)
                return Rcpp::wrap(1);
            else if (n < 5)
                return Rcpp::wrap(2);
            else if (n < 7)
                return Rcpp::wrap(3);
            else
                return Rcpp::wrap(4);
        }
        
        return Rcpp::wrap(static_cast<int>(PrimeCounting::PiPrime(n)));
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    int64_t result = PrimeCounting::MasterPrimeCount(n, nThreads, maxThreads);
    
    if (result > std::numeric_limits<int>::max())
        return Rcpp::wrap(static_cast<double>(result));
    else
        return Rcpp::wrap(static_cast<int>(result));
}

template <typename typeInt, typename typeReturn, typename typeRcpp>
SEXP GlueMotley(typeInt myMin, typeReturn myMax, bool isEuler,
                typeRcpp temp, bool keepNames, int nThreads, int maxThreads) {
    
    std::size_t myRange = (myMax - myMin) + 1;
    std::vector<typeReturn> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = myMin;
        for (std::size_t k = 0; retM <= myMax; ++retM, ++k)
            myNames[k] = retM;
    }
    
    if (isEuler) {
        std::vector<std::vector<typeInt>> tempList;
        typeRcpp EulerPhis(myRange);
        std::vector<typeInt> numSeq(myRange);
        MotleyPrimes::MotleyMaster(myMin, myMax, isEuler, EulerPhis,
                                   numSeq, tempList, nThreads, maxThreads);
        if (keepNames)
            EulerPhis.attr("names") = myNames;
        
        return EulerPhis;
    } else {
        std::vector<std::vector<typeInt>> 
            primeList(myRange, std::vector<typeInt>());
        typeRcpp tempRcpp;
        std::vector<typeInt> tempVec;
        MotleyPrimes::MotleyMaster(myMin, myMax, isEuler, tempRcpp,
                                   tempVec, primeList, nThreads, maxThreads);
        
        Rcpp::List myList = Rcpp::wrap(primeList);
        if (keepNames)
            myList.attr("names") = myNames;
        
        return myList;
    }
}

// [[Rcpp::export]]
SEXP MotleyContainer(SEXP Rb1, SEXP Rb2, bool isEuler, 
                     SEXP RNamed, SEXP RNumThreads, int maxThreads) {
    
    double bound1, bound2, myMin, myMax;
    const std::string namedObject = (isEuler) ? "namedVector" : "namedList";
    bool isNamed = CleanConvert::convertLogical(RNamed, namedObject);
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1");
    
    if (bound1 <= 0 || bound1 > Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2");
    }
    
    if (bound2 <= 0 || bound2 > Significand53)
        Rcpp::stop("bound2 must be a positive number less than 2^53");
    
    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }
    
    if (myMax < 2) {
        if (isEuler) {
            Rcpp::IntegerVector z(1, 1);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        } else {
            std::vector<std::vector<int>> trivialRet(1, std::vector<int>());
            Rcpp::List z = Rcpp::wrap(trivialRet);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        }
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    if (myMax > std::numeric_limits<int>::max()) {
        int64_t intMin = static_cast<int64_t>(myMin);
        Rcpp::NumericVector temp;
        return GlueMotley(intMin, myMax, isEuler, temp, isNamed, nThreads, maxThreads);
    } else {
        int32_t intMin = static_cast<int32_t>(myMin);
        int32_t intMax = static_cast<int32_t>(myMax);
        Rcpp::IntegerVector temp;
        return GlueMotley(intMin, intMax, isEuler, temp, isNamed, nThreads, maxThreads);
    }
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads, int maxCores, int maxThreads) {
    
    double bound1, bound2;
    int_fast64_t myMax, myMin;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1", false, false);
    
    if (bound1 <= 0 || bound1 > Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2", false, false);
    }
    
    if (bound2 <= 0 || bound2 > Significand53)
        Rcpp::stop("bound2 must be a positive number less than 2^53");
    
    if (bound1 > bound2) {
        myMax = static_cast<int_fast64_t>(std::floor(bound1));
        myMin = static_cast<int_fast64_t>(std::ceil(bound2));
    } else {
        myMax = static_cast<int_fast64_t>(std::floor(bound2));
        myMin = static_cast<int_fast64_t>(std::ceil(bound1));
    }
    
    if (myMax <= 1)
        return Rcpp::IntegerVector();
    
    if (myMin <= 2) myMin = 1;
    if (myMin == myMax) {++myMax;}
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    std::size_t numPrimes = 0u;
    std::vector<unsigned long int> runningCount;
    runningCount.push_back(0u);
    unsigned long int numSects = nThreads;
    bool Parallel = false;
    
    if (myMax > std::numeric_limits<int>::max()) {
        std::vector<std::vector<double>> primeList(numSects, std::vector<double>());
        std::vector<double> tempPrime;
        
        PrimeSieve::PrimeSieveMaster(myMin, myMax, tempPrime, primeList,
                                     Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::NumericVector primes(numPrimes);
            Rcpp::NumericVector::iterator priBeg = primes.begin();
            
            for (std::size_t i = 0; i < numSects; ++i)
                std::move(primeList[i].cbegin(), primeList[i].cend(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    } else {
        std::vector<std::vector<int_fast32_t>> primeList(numSects, std::vector<int_fast32_t>());
        std::vector<int_fast32_t> tempPrime;
        
        PrimeSieve::PrimeSieveMaster(myMin, myMax, tempPrime, primeList,
                                     Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::IntegerVector primes(numPrimes);
            Rcpp::IntegerVector::iterator priBeg = primes.begin();
            
            for (std::size_t i = 0; i < numSects; ++i)
                std::move(primeList[i].cbegin(), primeList[i].cend(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    }
}

