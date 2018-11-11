#include <Rcpp.h>
#include <cmath>
#include <array>
#include <libdivide.h>
#include <PrimesSegSieve.h>
#include <PhiTinyLookup.h>
#include <GetFacsUtils.h>
#include <CleanConvert.h>
#include <thread>

// "PrimeSieve" implements a simple segmented version of the Sieve of 
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve
// "MasterPrimeCount" is a simple variation of the highly optimized "pi_legendre.cpp"
// algorithm by Kim Walisch, which calculates the numbers of primes less than n using
// Legendre's formula. Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

// This is the largest multiple of 2*3*5*7 = 210
// that is less than 2^15 = 32768 = 32KB. This
// is the typical size of most CPU's L1 cache
const int L1CacheSize = 32760;
const unsigned long int wheelSize = 48;
const int maxVal210 = 210;
const int_fast64_t segmentSize = (int_fast64_t) L1CacheSize;
constexpr unsigned long int numSegments = L1CacheSize / maxVal210;

static const int_fast64_t wheel210[wheelSize] = {10, 2, 4, 2, 4, 6, 2, 6, 4, 2,
                                                  4, 6, 6, 2, 6, 4, 2, 6, 4, 6,
                                                  8, 4, 2, 4, 2, 4, 8, 6, 4, 6,
                                                  2, 4, 6, 2, 6, 6, 4, 2, 4, 6,
                                                  2, 6, 4, 2, 4, 2, 10, 2};

const unsigned long int lastSmlPri = smallPrimeBase[(sizeof(smallPrimeBase) / sizeof(smallPrimeBase[0])) - 1];

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

const std::vector<double> percInc = {0.25, 0.116, 0.103, 0.0850, 0.0712,
                                     0.0614, 0.0538, 0.0480, 0.0431, 
                                     0.0392, 0.0360, 0.0332, 0.0309};

const std::vector<double> cutPoints = {40000.0, 120000.0, 1000000.0, 10000000.0,
                                       100000000.0, 1000000000.0, 10000000000.0,
                                       100000000000.0, 1000000000000.0,
                                       10000000000000.0, 100000000000000.0,
                                       1000000000000000.0, 10000000000000000.0};

// The following function is based off of the prime number theorem
std::size_t EstimatePrimeCount(double minNum, double maxNum) {
    std::vector<double>::const_iterator it = std::upper_bound(cutPoints.begin(),
                                                              cutPoints.end(),
                                                              maxNum);
    std::size_t myIndex = it - cutPoints.begin();
    double dblRes = std::ceil((maxNum / log(maxNum)) * (1 + percInc[myIndex]));
    
    if (minNum > 1000)
        dblRes -= std::floor((minNum / log(minNum)) * (1 + percInc[myIndex]));
    
    std::size_t result = dblRes;
    return result;
}

template <typename typeReturn, typename typePrime>
void PrimeSieve(typePrime minNum, typePrime maxNum,
                const std::vector<int_fast64_t> &smallPrimes,
                int sqrtBound, std::vector<typeReturn> &myPrimes) {
    
    int_fast64_t segSize = segmentSize;
    unsigned long int numSegs = numSegments;
    std::size_t myReserve = EstimatePrimeCount((double) minNum, (double) maxNum);
    myPrimes.reserve(myReserve);
    
    if (maxNum <= lastSmlPri) {
        std::size_t ind = 0;
        for (; smallPrimeBase[ind] < minNum; ++ind) {}
        
        for (; smallPrimeBase[ind] <= maxNum; ++ind) 
            myPrimes.push_back((typeReturn) smallPrimeBase[ind]);
    } else {
        if (minNum < 13) {
            std::size_t ind = 0;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}
            
            for (; smallPrimeBase[ind] < 10; ++ind)
                myPrimes.push_back((typeReturn) smallPrimeBase[ind]);
        }
        
        // Ensure segSize is greater than sqrt(n)
        // as well as multiple of L1CacheSize
        if (segSize < sqrtBound) {
            numSegs = (unsigned long int) ceil((double) sqrtBound / maxVal210);
            segSize = numSegs * maxVal210;
        }
        
        std::vector<int_fast64_t> nextStrt;
        int_fast64_t flrMaxNum = segSize * floor((double) maxNum / segSize);
        
        // vector used for sieving
        std::vector<char> sieve(segSize, 1);
        sieve[1] = 0;
        
        int_fast64_t sqrPrime = (int_fast64_t) (3 * 3);
        int_fast64_t lowerBnd = 0, upperBnd = segSize, myNum = 1;
        int p = 1;
        
        if (minNum > 2) {
            lowerBnd = (typePrime) (segSize * floor((double) minNum / segSize));
            upperBnd = std::min(lowerBnd + segSize, (int_fast64_t) maxNum);
            myNum += lowerBnd;
            int_fast64_t myStart;
            
            for (; sqrPrime <= upperBnd; ++p) {
                if (lowerBnd > sqrPrime) {
                    int_fast64_t remTest = lowerBnd % smallPrimes[p - 1];
                    if (remTest == 0) {
                        myStart = smallPrimes[p - 1];
                    } else {
                        myStart = smallPrimes[p - 1] - remTest;
                        if ((myStart % 2) == 0) {myStart += smallPrimes[p - 1];}
                    }
                } else {
                    myStart = sqrPrime - lowerBnd;
                }
                
                nextStrt.push_back(myStart);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            if (upperBnd < flrMaxNum) {
                for (std::size_t q = 0; q < numSegs; ++q) {
                    for (std::size_t w = 0; w < wheelSize; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((typeReturn) myNum);
                        myNum += wheel210[w];
                    }
                }
            } else {
                for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                    for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((typeReturn) myNum);
                        myNum += wheel210[w];
                    }
                }
            }
            
            std::fill(sieve.begin(), sieve.end(), 1);
            lowerBnd += segSize;
        }
        
        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            // current segment = interval [lowerBnd, upperBnd]
            upperBnd = lowerBnd + segSize;
            
            // sieve the current segment && sieving primes <= sqrt(upperBnd)
            for (; sqrPrime <= upperBnd; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < numSegs; ++q) {
                for (std::size_t w = 0; w < wheelSize; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((typeReturn) myNum);
                    
                    myNum += wheel210[w];
                }
            }
            
            std::fill(sieve.begin(), sieve.end(), 1);
        }
        
        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
            }
            
            for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((typeReturn) myNum);
                    
                    myNum += wheel210[w];
                }
            }
        }
    }
}

template <typename typePrime>
std::vector<typePrime> sqrtBasePrimes(int sqrtBound, bool bAddZero,
                                      bool bAddExtraPrime, bool bAddTwo) {
    
    std::vector<typePrime> smallPrimes;
    
    if (sqrtBound < lastSmlPri) {
        if (bAddZero) smallPrimes.push_back(0);
        unsigned long int ind = (bAddTwo) ? 0 : 1;
        
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            smallPrimes.push_back((typePrime) smallPrimeBase[ind]);
        
        if (bAddExtraPrime)
            smallPrimes.push_back((typePrime) smallPrimeBase[ind]);
    } else {
        int sqrtSqrtBound = (int) std::sqrt(sqrtBound);
        std::vector<int_fast64_t> smallSmlPrimes;
        unsigned long int ind = 1;
        
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            smallSmlPrimes.push_back((typePrime) smallPrimeBase[ind]);
        
        smallSmlPrimes.push_back((typePrime) smallPrimeBase[ind]);
        
        // The number, 225, comes from the observation that the largest prime
        // gap less than 100 million is 219 @ 47,326,693. This is important
        // because we need to guarantee that we obtain the smallest prime
        // greater than sqrt(n). We know that the largest number that this
        // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266
        
        typePrime myLower = 3, myUpper = sqrtBound;
        if (bAddExtraPrime) {myUpper += 225;}
        
        std::size_t sqrtReserve = EstimatePrimeCount(1.0, (double) myUpper);
        smallPrimes.reserve(sqrtReserve);
        
        if (bAddZero) smallPrimes.push_back(0);
        if (bAddTwo) {myLower = 1;}
        
        PrimeSieve(myLower, myUpper, smallSmlPrimes, sqrtSqrtBound, smallPrimes);
    }
    
    return smallPrimes;
}

// PiPrime is very similar to the AllPrimesC only we are not
// considering a range. That is, we are only concerned with finding
// the number of primes less than maxNum. We are also only counting
// primes instead of generating them.
int64_t PiPrime (int64_t maxNum) {
    
    int_fast64_t segSize = segmentSize;
    unsigned long int numSegs = numSegments;
    int sqrtBound = (int) std::sqrt((double) maxNum);
    
    // the wheel already has the first 4 primes marked as
    // false, so we need to account for them here. N.B.
    // the calling functions has checks for cases where
    // maxNum is less than 11, so need to check here.
    int64_t count = 4;
    
    if (segSize < sqrtBound) {
        numSegs = (unsigned long int) ceil((double) sqrtBound / maxVal210);
        segSize = numSegs * maxVal210;
    }
    
    std::vector<int64_t> smallPrimes, nextStrt;
    int64_t flrMaxNum = segSize * floor((double) maxNum / segSize);
    
    int i = 1;
    for (; smallPrimeBase[i] <= sqrtBound; ++i)
        smallPrimes.push_back(smallPrimeBase[i]);
    
    smallPrimes.push_back(smallPrimeBase[i]);
    std::vector<char> sieve(segSize, 1);
    sieve[1] = 0;
    
    int64_t sqrPrime = (int64_t) (3 * 3);
    int64_t lowerBnd = 0, upperBnd = segSize, myNum = 1;
    int p = 1;
    
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
        
        for (std::size_t q = 0; q < numSegs; ++q) {
            for (std::size_t w = 0; w < wheelSize; ++w) {
                if (sieve[myNum - lowerBnd])
                    ++count;
                myNum += wheel210[w];
            }
        }
        
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
        
        for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
            for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                if (sieve[myNum - lowerBnd])
                    ++count;
                myNum += wheel210[w];
            }
        }
    }
    
    return count;
}

const int MAX_A = 100;
const int phiTinySize = 6;
std::array<std::vector<uint16_t>, MAX_A> phiCache;

// Increment MAX_A, so we can have easier access
// to indexing.  E.g when a = 3 (i.e. the number)
// of primes is 3), instead of accessing the third
// entry in phiTiny like phiTiny[a - 1], we simply
// use "a" as our index (i.e. phiTiny[a]). The
// first entry in phiTiny will be empty. The same
// goes for phiPrimes and phiPi.
std::array<std::vector<int16_t>, phiTinySize + 1> phiTiny;
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
        if (x >= phiCache[a].size())
            phiCache[a].resize(x + 1, 0);
        
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

//[[Rcpp::export]]
SEXP MasterPrimeCount (SEXP Rn) {
    double dblNum;
    
    switch(TYPEOF(Rn)) {
        case REALSXP: {
            dblNum = Rcpp::as<double>(Rn);
            break;
        }
        case INTSXP: {
            dblNum = Rcpp::as<double>(Rn);
            break;
        }
        default: {
            Rcpp::stop("n must be of type numeric or integer");
        }
    }

    if (dblNum < 1 || dblNum > Significand53)
        Rcpp::stop("n must be a positive number less than 2^53");

    int64_t n = (int64_t) dblNum;
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
        
        return Rcpp::wrap((int) PiPrime(n));
    }
    
    // populate phiTiny... phiTiny[0] should not be accessed
    phiTiny[1].resize(1 + 1);
    phiTiny[1][0] = 0;
    phiTiny[1][1] = 1;
    int16_t arr6[] = {0,1,1,1,1,2};
    std::vector<int16_t> phi6 (arr6, arr6 + sizeof(arr6) / sizeof(arr6[0]) );
    phiTiny[2] = phi6;
    std::vector<int16_t> phi30 (arr30, arr30 + sizeof(arr30) / sizeof(arr30[0]) );
    phiTiny[3] = phi30;
    std::vector<int16_t> phi210 (arr210, arr210 + sizeof(arr210) / sizeof(arr210[0]) );
    phiTiny[4] = phi210;
    std::vector<int16_t> phi2310 (arr2310, arr2310 + sizeof(arr2310) / sizeof(arr2310[0]) );
    phiTiny[5] = phi2310;
    
    std::vector<int16_t> phi30030;
    phi30030.reserve(30031);
    phi30030.push_back(0);
    
    for (std::size_t i = 1; i <= 5760; ++i)
        for (int16_t j = 0; j < arr30030freq[i]; ++j)
            phi30030.push_back(i);
    
    phiTiny[6] = phi30030;
    
    int64_t sqrtBound = (int64_t) std::sqrt((double) n);
    phiPrimes = sqrtBasePrimes<int_fast64_t>(sqrtBound, true, false, true);
    
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
    
    if (int64result > INT_MAX) {
        double dblResult = (double) int64result;
        return Rcpp::wrap(dblResult);
    }
    
    int intResult = (int) int64result;
    return Rcpp::wrap(intResult);
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
        std::vector<typeInt> primes = sqrtBasePrimes<typeInt>(sqrtBound, false, false, true);
        typename std::vector<typeInt>::iterator p, primesEnd;
        std::vector<int> myMemory(myRange, 1);
        std::vector<int>::iterator myMalloc;
        primesEnd = primes.end();
        
        for (p = primes.begin(); p < primesEnd; ++p) {
            limit = (unsigned long int) trunc(myLogN/log((double)*p));
            if (m < 2) {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
                    for (j = (myStep - 1); j < n; j += myStep)
                        ++myMemory[j];
                }
            } else {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
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
                limit = (unsigned long int) trunc(myLogN/log((double)*p));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
                    
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
                limit = (unsigned long int) trunc(myLogN/log((double)*p));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double)*p, (double) i);
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
    
    typeInt j, priTypeInt, myNum = m;
    double myLogN = log((double) n);
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
    int sqrtBound = floor(sqrt((double) n));
    typename std::vector<typeInt>::iterator p;
    
    if (m < 2) {
        std::vector<int_fast64_t> smallPrimes = sqrtBasePrimes<int_fast64_t>(sqrtBound, false, true, false);
        PrimeSieve((typeInt) 1, n, smallPrimes, sqrtBound, primes);
        
        for (p = primes.begin(); p < primes.end(); ++p) {
            libdivide::divider<typeInt> fastDiv(*p);
            for (j = (*p - 1); j < n; j += *p) {
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
        }
    } else if (n > 3) {
        primes = sqrtBasePrimes<typeInt>(sqrtBound, false, false, true);
        typeInt myStart, myStep;
    
        for (p = primes.begin(); p < primes.end(); ++p) {
            limit = (unsigned long int) trunc(myLogN / log((double) *p));
            priTypeInt = *p;
            myStart = getStartIndexPowP(m, *p, priTypeInt);
            libdivide::divider<typeInt> fastDiv(priTypeInt);
    
            for (j = myStart; j < myRange; j += *p) {
                numSeq[j] /= fastDiv;
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
    
            for (std::size_t i = 2; i <= limit; ++i) {
                myStep = (typeInt) pow((double)*p, i);
                myStart = getStartIndexPowP(m, myStep, priTypeInt);
                
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
unsigned int TotalNumThreads() {
    return std::thread::hardware_concurrency();
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2, SEXP RIsList, SEXP RIsEuler,
                       SEXP RNamed, SEXP RNumThreads) {
    double bound1, bound2, myMax, myMin;
    bool Parallel = false, isList = false, isEuler = false, isNamed = false;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1 must be of type numeric or integer", false);
    
    isList = Rcpp::as<bool>(RIsList);
    isEuler = Rcpp::as<bool>(RIsEuler);
    isNamed = Rcpp::as<bool>(RNamed);
    
    if (bound1 <= 0 || bound1 > Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2 must be of type numeric or integer", false);
    }
    
    if (bound2 <= 0 || bound2 > Significand53)
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
    } else if (isEuler) {
        if (myMax <= 1) {
            Rcpp::IntegerVector z(1, 1);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        }
        
        if (myMax >= INT_MAX)
            return EulerPhiSieveCpp<Rcpp::NumericVector>((int64_t) myMin, (double) myMax, isNamed);
        
        return EulerPhiSieveCpp<Rcpp::IntegerVector>((int32_t) myMin, (int32_t) myMax, isNamed);
    } else {
        if (myMax <= 1)
            return Rcpp::IntegerVector();
        
        if (myMin <= 2) myMin = 1;
        if (myMin == myMax) {++myMax;}
        
        int numThreads, sqrtBound = (int) std::sqrt(myMax);
        std::vector<int_fast64_t> smallPrimes = sqrtBasePrimes<int_fast64_t>(sqrtBound, false, true, false);
        int_fast64_t myRange = (int_fast64_t) (myMax - myMin);
        int totalThreads = (int) TotalNumThreads();
        
        if ((myRange < 200000) || (totalThreads < 2)) {
            Parallel = false;
        } else if (!Rf_isNull(RNumThreads)) {
            Parallel = true;
            CleanConvert::convertPrimitive(RNumThreads, numThreads, "nThreads must be of type numeric or integer");
            if (numThreads > totalThreads)
                numThreads = totalThreads;
            if (numThreads < 2)
                Parallel = false;
        }
        
        if (Parallel) {
            std::vector<std::thread> myThreads;
            std::size_t numPrimes = 0, count = 0;
            std::vector<unsigned long int> sectionSize(numThreads);
            
            if (myMax >= INT_MAX) {
                std::vector<std::vector<double> > primeList(numThreads, std::vector<double>());
                int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
                int_fast64_t upperBnd = myMin + stepSize - 1;

                for (int j = 0; j < (numThreads - 1); ++j) {
                    myThreads.emplace_back(PrimeSieve<double, int_fast64_t>, lowerBnd,
                                           upperBnd, std::ref(smallPrimes), sqrtBound, std::ref(primeList[j]));
                    lowerBnd += stepSize;
                    upperBnd += stepSize;
                }

                myThreads.emplace_back(PrimeSieve<double, int_fast64_t>, lowerBnd,
                                       (int_fast64_t) myMax, std::ref(smallPrimes),
                                       sqrtBound, std::ref(primeList[numThreads - 1]));

                for (auto& thr: myThreads)
                    thr.join();

                for (int i = 0; i < numThreads; ++i) {
                    numPrimes += primeList[i].size();
                    sectionSize[i] = primeList[i].size();
                }

                if (numPrimes < INT_MAX) {
                    Rcpp::NumericVector primes(numPrimes);

                    for (int i = 0; i < numThreads; ++i)
                        for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                            primes[count] = primeList[i][j];

                    return primes;
                } else {
                    return Rcpp::wrap(primeList);
                }
            }
            
            std::vector<std::vector<int_fast32_t> > primeList(numThreads, std::vector<int_fast32_t>());
            int_fast32_t lowerBnd = myMin, stepSize = myRange / numThreads;
            int_fast32_t upperBnd = myMin + stepSize - 1;
            
            for (int j = 0; j < (numThreads - 1); ++j) {
                myThreads.emplace_back(PrimeSieve<int_fast32_t, int_fast32_t>, lowerBnd,
                                       upperBnd, std::ref(smallPrimes), sqrtBound, std::ref(primeList[j]));
                lowerBnd += stepSize;
                upperBnd += stepSize;
            }
            
            myThreads.emplace_back(PrimeSieve<int_fast32_t, int_fast32_t>, lowerBnd,
                                   (int_fast32_t) myMax, std::ref(smallPrimes),
                                   sqrtBound, std::ref(primeList[numThreads - 1]));

            for (auto& thr: myThreads)
                thr.join();
            
            for (int i = 0; i < numThreads; ++i) {
                numPrimes += primeList[i].size();
                sectionSize[i] = primeList[i].size();
            }
            
            Rcpp::IntegerVector primes(numPrimes);

            for (int i = 0; i < numThreads; ++i)
                for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                    primes[count] = primeList[i][j];

            return primes;
        }
        
        if (myMax >= INT_MAX) {
            std::vector<double> primes;
            PrimeSieve((int_fast64_t) myMin, (int_fast64_t) myMax, smallPrimes, sqrtBound, primes);
            return Rcpp::wrap(primes);
        }
        
        std::vector<int_fast32_t> primes;
        PrimeSieve((int_fast32_t) myMin, (int_fast32_t) myMax, smallPrimes, sqrtBound, primes);
        return Rcpp::wrap(primes);
    }
}
