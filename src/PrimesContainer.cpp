#include "EratosthenesSieve.h"
#include "CleanConvert.h"
#include "PhiTinyLookup.h"
#include <libdivide.h>
#include <array>

namespace PrimeCounting {
    
    // This is the largest multiple of 2*3*5*7 = 210
    // that is less than 2^15 = 32768 = 32KB. This
    // is the typical size of most CPU's L1 cache
    constexpr int_fast64_t AlmostL1Cache = 32760;

    constexpr unsigned long int SZ_WHEEL210 = 48;
    constexpr unsigned long int NUM210 = 210;
    constexpr std::size_t N_WHEELS_PER_SEG = static_cast<std::size_t>(AlmostL1Cache / NUM210);
    
    static const int_fast64_t ARR_WHEEL210[SZ_WHEEL210] = {
        10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
        2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2};
    
    // PiPrime is very similar to the PrimeSieveSmall only we are not
    // considering a range. That is, we are only concerned with finding
    // the number of primes less than maxNum. We are also only counting
    // primes instead of generating them.
    int64_t PiPrime (int64_t maxNum) {
        
        int_fast64_t segSize = AlmostL1Cache;
        unsigned long int numSegs = N_WHEELS_PER_SEG;
        std::size_t szWheel210 = SZ_WHEEL210;
        int sqrtBound = static_cast<int>(std::sqrt(maxNum));
        
        // the wheel already has the first 4 primes marked as
        // false, so we need to account for them here. N.B.
        // the calling functions has checks for cases where
        // maxNum is less than 11, so no need to check here.
        int64_t count = 4;
        
        std::vector<int64_t> smallPrimes, nextStrt;
        int64_t flrMaxNum = segSize * std::floor(maxNum / segSize);
        
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
            
            for (std::size_t q = 0; q < numSegs; ++q)
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
            
            for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += ARR_WHEEL210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        ++count;
        }
        
        return count;
    }

    const int MAX_A = 100;
    std::array<std::vector<uint16_t>, MAX_A> phiCache;
    
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
    const std::array<int, 7> myTinyPrimes = {{0, 2, 3, 5, 7, 11, 13}};
    const std::array<int, 7> primeProds = {{1, 2, 6, 30, 210, 2310, 30030}};
    const std::array<int, 7> myTotients = {{1, 1, 2, 8, 48, 480, 5760}};
    const std::array<int, 13> myTinyPi = {{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5}};
    
    inline int64_t phiTinyCalc (int64_t x, int64_t a) {
        int64_t pp = primeProds[a];
        return (x / pp) * myTotients[a] + phiTiny[a][x % pp];
    }
    
    void updateCache (uint64_t x, uint64_t a, int64_t mySum) {
        if (a < phiCache.size() &&
            x <= std::numeric_limits<uint16_t>::max()) {
            if (x >= phiCache[a].size()) phiCache[a].resize(x + 1, 0);
            phiCache[a][x] = (uint16_t) std::abs(mySum);
        }
    }
    
    inline bool isCached (uint64_t x, uint64_t a) {
        return a < phiCache.size() &&
               x < phiCache[a].size() &&
               phiCache[a][x];
    }
    
    inline bool isPix(int64_t x, int64_t a) {
        return x < (int) phiPi.size() &&
               x < (phiPrimes[a + 1] * phiPrimes[a + 1]);
    }
    
    inline int64_t getStrt(int64_t y) {
        if (y >= myTinyPrimes.back())
            return phiTinySize;
        else
            return myTinyPi[y];
    }
    
    template <int SIGN>
    int64_t phiWorker (int64_t x, int64_t a) {
        if (x <= phiPrimes[a])
            return SIGN;
        else if (a <= phiTinySize)
            return phiTinyCalc(x, a) * SIGN;
        else if (isPix(x, a))
            return (phiPi[x] - a + 1) * SIGN;
        else if (isCached(x, a))
            return phiCache[a][x] * SIGN;
            
        int64_t sqrtx = (int64_t) std::sqrt((double) x);
        int64_t piSqrtx = a;
        int64_t strt = getStrt(sqrtx);
        int64_t mySum = 0;
        
        if (sqrtx < (int) phiPi.size())
            piSqrtx = std::min((int64_t) phiPi[sqrtx], a);
        
        mySum += (piSqrtx - a) * SIGN;
        mySum += phiTinyCalc(x, strt) * SIGN;
        
        for (int64_t i = strt; i < piSqrtx; ++i) {
            int64_t x2 = x / phiPrimes[i + 1];
            
            if (isPix(x2, i))
                mySum += (phiPi[x2] - i + 1) * -SIGN;
            else
                mySum += phiWorker<-SIGN>(x2, i);
        }
        
        updateCache(x, a, mySum);
        return mySum;
    }
    
    int64_t phiMaster (int64_t x, int64_t a) {
        
        int64_t sqrtx = (int64_t) std::sqrt((double) x);
        int64_t piSqrtx = std::min((int64_t) phiPi[sqrtx], a);
        int64_t strt = getStrt(sqrtx);
        int64_t mySum = phiTinyCalc(x, strt) + piSqrtx - a;
        
        for (int64_t i = strt; i < piSqrtx; ++i)
            mySum += phiWorker<-1>(x / phiPrimes[i + 1], i);
        
        return mySum;
    }
    
    // All values verified by Kim Walisch's primecount library
    // 10^9  -->> 50,847,534
    // 10^10 -->> 455,052,511
    // 10^11 -->> 4,118,054,813
    // 10^12 -->> 37,607,912,018
    // 10^13 -->> 346,065,536,839
    // 10^14 -->> 3,204,941,750,802
    // 10^15 -->> 29,844,570,422,669
    // MAX VALUE (2^53 - 1) -->> 252,252,704,148,404
    
    int64_t MasterPrimeCount (int64_t n) {
        
        int64_t sqrtBound = static_cast<int64_t>(std::sqrt(n));
        std::vector<int64_t> resetPhiPrimes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, true, false, true, resetPhiPrimes);
        phiPrimes = resetPhiPrimes;
        
        phiPi.resize(sqrtBound + 1);
        int64_t count = 0;
        int64_t maxPrime = phiPrimes.back();
        
        for (int64_t i = 1; i <= maxPrime; ++i) {
            if (i >= phiPrimes[count + 1])
                ++count;
            phiPi[i] = count;
        }
        
        for (int64_t i = (maxPrime + 1); i <= sqrtBound; ++i)
            phiPi[i] = count;
        
        int64_t piSqrt = PiPrime(sqrtBound);
        int64_t phiSqrt = phiMaster(n, piSqrt);
        int64_t int64result = piSqrt + phiSqrt - 1;
        
        return int64result;
    }
}

//[[Rcpp::export]]
SEXP PrimeCountRcpp (SEXP Rn) {
    double dblNum;
    CleanConvert::convertPrimitive(Rn, dblNum, "n must be of type numeric or integer", false);

    if (dblNum < 1 || dblNum > PrimeSieve::Significand53)
        Rcpp::stop("n must be a positive number less than 2^53");

    int64_t n = dblNum;
    
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
    
    int64_t result = PrimeCounting::MasterPrimeCount(n);
    
    if (result > INT_MAX)
        return Rcpp::wrap(static_cast<double>(result));
    else
        return Rcpp::wrap(static_cast<int>(result));
}

// This function is slightly different than the getStartingIndex
// in the DivisorsContainer.cpp file. The step passed in this
// function is a power of prime and requires an additional check
// (i.e. else if (myPrime < lowerB)).
template <typename typeInt>
inline typeInt getStartIndexPowP(typeInt lowerB, typeInt step, 
                                 typeInt myPrime) {
    
    typeInt retStrt, remTest = lowerB % step;
    
    if (remTest == 0) {
        retStrt = 0;
    } else if (myPrime < lowerB) {
        retStrt = step - remTest;
    } else {
        retStrt = step - lowerB;
    }
    
    return retStrt;
}

template <typename typeInt, typename typeReturn>
Rcpp::List PrimeFactorizationSieve(typeInt m, typeReturn retN, bool keepNames) {
    
    typeInt n = (typeInt) retN;
    typeInt myRange = n;
    myRange += (1 - m);
    
    std::vector<std::vector<typeReturn> > 
        MyPrimeList(myRange, std::vector<typeReturn>());
    
    typename std::vector<std::vector<typeReturn> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    typeInt sqrtBound = (typeInt) floor(sqrt((double)n));
    
    typeInt j, myStep, myStart, myNum = m;
    double myLogN = log((double)n);
    unsigned long int limit;
    
    std::vector<typeReturn> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = (typeReturn) m;
        for (std::size_t k = 0; retM <= retN; ++retM, ++k)
            myNames[k] = retM;
    }
    
    if (n > 3) {
        std::vector<typeInt> primes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, false, false, true, primes);
        typename std::vector<typeInt>::iterator p, primesEnd;
        std::vector<int> myMemory(myRange, 1);
        std::vector<int>::iterator myMalloc;
        primesEnd = primes.end();
        
        for (p = primes.begin(); p < primesEnd; ++p) {
            limit = static_cast<unsigned long int>(trunc(myLogN /
                                    log(static_cast<double>(*p))));
            if (m < 2) {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = static_cast<typeInt>(pow(*p, i));
                    for (j = (myStep - 1); j < n; j += myStep)
                        ++myMemory[j];
                }
            } else {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = static_cast<typeInt>(pow(*p, i));
                    myStart = getStartIndexPowP(m, myStep, *p);
                    for (j = myStart; j < myRange; j += myStep)
                        ++myMemory[j];
                }
            }
        }
            
        if (myNum < 2){
            ++myNum;
            it2d = MyPrimeList.begin() + 1;
            myMalloc = myMemory.begin() + 1;
        } else {
            it2d = MyPrimeList.begin();
            myMalloc = myMemory.begin();
        }
        
        for (; it2d < itEnd; ++it2d, ++myNum, ++myMalloc) {
            it2d->reserve(*myMalloc);
            it2d->push_back((typeReturn) myNum);
        }
    
        if (m < 2) {
            for (p = primes.begin(); p < primesEnd; ++p) {
                limit = static_cast<unsigned long int>(trunc(myLogN /
                                    log(static_cast<double>(*p))));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = static_cast<typeInt>(pow(*p, i));
                    
                    for (j = (myStep - 1); j < n; j += myStep) {
                        if (MyPrimeList[j].back() > *p) {
                            myNum = (typeInt) MyPrimeList[j].back();
                            myNum /= fastDiv;
                            MyPrimeList[j].back() = (typeReturn) myNum;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, (typeReturn) *p);
                        }
                    }
                }
            }
        } else {
            for (p = primes.begin(); p < primesEnd; ++p) {
                limit = static_cast<unsigned long int>(trunc(myLogN /
                                        log(static_cast<double>(*p))));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = static_cast<typeInt>(pow(*p, i));
                    myStart = getStartIndexPowP(m, myStep, *p);
    
                    for (j = myStart; j < myRange; j += myStep) {
                        if (MyPrimeList[j].back() > *p) {
                            myNum = (typeInt) MyPrimeList[j].back();
                            myNum /= fastDiv;
                            MyPrimeList[j].back() = (typeReturn) myNum;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, (typeReturn) *p);
                        }
                    }
                }
            }
        }
    } else { // edge case where m,n = 2 or 3
        int strt = 0;
        myNum = m;
        if (m == 1) {
            ++strt;
            ++myNum;
        }
        for (int i = strt; i < myRange; ++i, ++myNum)
            MyPrimeList[i].push_back(myNum);
    }
    
    Rcpp::List myList = Rcpp::wrap(MyPrimeList);
    if (keepNames)
        myList.attr("names") = myNames;
    
    return myList;
}

template <typename typeRcpp, typename typeInt, typename typeReturn>
typeRcpp EulerPhiSieveCpp (typeInt m, typeReturn retN, bool keepNames) {
    
    typeInt n = (typeInt) retN;
    typeInt myRange = n;
    myRange += (1 - m);
    
    typeInt j, myNum = m;
    double myLogN = log(static_cast<double>(n));
    unsigned long int limit;
    std::vector<typeReturn> EulerPhis(myRange);
    std::vector<typeInt> numSeq(myRange);
    
    std::vector<typeReturn> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = (typeReturn) m;
        for (std::size_t k = 0; retM <= retN; ++retM, ++k)
            myNames[k] = retM;
    }
    
    typeReturn retM = (typeReturn) m;
    for (std::size_t i = 0; retM <= retN; ++retM, ++i) {
        EulerPhis[i] = retM;
        numSeq[i] = (typeInt) retM;
    }
    
    std::vector<typeInt> primes;
    int sqrtBound = static_cast<int>(std::floor(sqrt(n)));
    typename std::vector<typeInt>::iterator p;
    
    if (m < 2) {
        bool Parallel = false;
        std::vector<std::vector<typeInt>> primeList;
        int_fast64_t one64 = 1, n64 = static_cast<int_fast64_t>(n);
        PrimeSieve::PrimeMaster(one64, n64, primes, primeList, Parallel);

        for (p = primes.begin(); p < primes.end(); ++p) {
            libdivide::divider<typeInt> fastDiv(*p);
            for (j = (*p - 1); j < n; j += *p) {
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
        }
    } else if (n > 3) {
        PrimeSieve::sqrtBigPrimes(sqrtBound, false, false, true, primes);
        typeInt myStart, myStep;
    
        for (p = primes.begin(); p < primes.end(); ++p) {
            limit = static_cast<unsigned long int>(std::trunc(myLogN /
                                    log(static_cast<double>(*p))));
            typeInt myP = static_cast<typeInt>(*p);
            myStart = getStartIndexPowP(m, myP, myP);
            libdivide::divider<typeInt> fastDiv(*p);
    
            for (j = myStart; j < myRange; j += myP) {
                numSeq[j] /= fastDiv;
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
    
            for (std::size_t i = 2; i <= limit; ++i) {
                myStep = (typeInt) pow(*p, i);
                myStart = getStartIndexPowP(m, myStep, myP);
                
                for (j = myStart; j < myRange; j += myStep)
                    numSeq[j] /= fastDiv;
            }
        }

        for (typeInt i = 0; i < myRange; ++i) {
            if (numSeq[i] > 1) {
                myNum = (typeInt) EulerPhis[i];
                myNum /= numSeq[i];
                EulerPhis[i] -= (typeReturn) myNum;
            }
        }
    } else { // edge case where m,n = 2 or 3
        for (int i = 0; i < myRange; ++i)
            --EulerPhis[i];
    }
    
    typeRcpp myVector = Rcpp::wrap(EulerPhis);
    if (keepNames)
        myVector.attr("names") = myNames;
        
    return myVector;
}

// [[Rcpp::export]]
SEXP MotleyPrimes(SEXP Rb1, SEXP Rb2, SEXP RIsList, 
                  SEXP RNamed, SEXP RNumThreads) {
    
    double bound1, bound2;
    int_fast64_t myMax, myMin;
    bool isList = false, isNamed = false;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1 must be of type numeric or integer", false);
    
    isList = Rcpp::as<bool>(RIsList);
    isNamed = Rcpp::as<bool>(RNamed);
    
    if (bound1 <= 0 || bound1 > PrimeSieve::Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2 must be of type numeric or integer", false);
    }
    
    if (bound2 <= 0 || bound2 > PrimeSieve::Significand53)
        Rcpp::stop("bound2 must be a positive number less than 2^53");
    
    if (bound1 > bound2) {
        myMax = bound1;
        myMin = bound2;
    } else {
        myMax = bound2;
        myMin = bound1;
    }
    
    myMin = ceil(myMin);
    myMax = floor(myMax);
        
    if (isList) {
        if (myMax < 2) {
            std::vector<std::vector<int> > trivialRet(1, std::vector<int>());
            Rcpp::List z = Rcpp::wrap(trivialRet);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        }
        
        if (myMax >= INT_MAX)
            return PrimeFactorizationSieve((int64_t) myMin, (double) myMax, isNamed);
        
        return PrimeFactorizationSieve((int32_t) myMin, (int32_t) myMax, isNamed);
    } else {
        if (myMax <= 1) {
            Rcpp::IntegerVector z(1, 1);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        }
        
        if (myMax >= INT_MAX)
            return EulerPhiSieveCpp<Rcpp::NumericVector>((int64_t) myMin, (double) myMax, isNamed);
        
        return EulerPhiSieveCpp<Rcpp::IntegerVector>((int32_t) myMin, (int32_t) myMax, isNamed);
    }
}
