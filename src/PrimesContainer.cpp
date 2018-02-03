#include <Rcpp.h>
#include <math.h>
#include <stdint.h>
#include "PrimesSegSieve.h"
#include "WheelInt30030.h"

using namespace Rcpp;

// Below are two functions for generating prime numbers quickly.
// "AllPrimesCpp" implements a segmented version of the Sieve of
// Eratosthenes (original implementation authored by Kim Walisch).
// An overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// The official github repo is:
//                      https://github.com/kimwalisch/primesieve
// "RangePrimeC" will return all primes within a given range.

const double Significand53 = 9007199254740991;

// 2*3*5*7*11*13...the primes used to generate the wheel
const int L1CacheSize = 30030;

// Number of elements in the static array wheelDiff30030
#define initialWheelSize (sizeof(wheelDiff30030) / sizeof(wheelDiff30030[0]))

template <typename typeRcpp, typename typePrime>
typeRcpp AllPrimesCpp (typePrime m, typePrime n) {

    int sqrtBound = (int) floor(sqrt((double)n));
    int segmentSize = L1CacheSize;
    int numSegments = 1;
    std::vector<typePrime> myPrimes;
    typePrime minNum = m, maxNum = n;
    typePrime myRange = maxNum - minNum;
    
    // Percentages obtained here: https://en.wikipedia.org/wiki/Prime-counting_function
    typePrime myReserve;
    if (maxNum < 100000) {
        myReserve = (typePrime) floor((double) 2*myRange/log((double)myRange));
    } else if (maxNum < 10000000) {
        myReserve = (typePrime) floor((double) 1.1*myRange/log((double)myRange));
    } else if (maxNum < 100000000) {
        myReserve = (typePrime) floor((double) 1.07*myRange/log((double)myRange));
    } else if (maxNum < 1000000000) {
        myReserve = (typePrime) floor((double) 1.06*myRange/log((double)myRange));
    } else {
        myReserve = (typePrime) floor((double) 1.053*myRange/log((double)myRange));
    }
    
    myPrimes.reserve(myReserve);
    
    uint_fast64_t sqrPrime, lastP, upperBnd;
    lastP = 2;
    std::size_t i = 0;
    
    if (maxNum < segmentSize || minNum <= 13)
        for (; lastP < minNum; i++)
            lastP += primesDiffSeg[i];
    
    if (maxNum < segmentSize) {
        for (; lastP <= maxNum; i++) {
            myPrimes.push_back(lastP);
            lastP += primesDiffSeg[i];
        }
        return wrap(myPrimes);
    } else if (minNum <= 13) {
        for (; lastP < 15; i++) {
            myPrimes.push_back(lastP);
            lastP += primesDiffSeg[i];
        }
    }
    
    // ensure segmentSize is greater than sqrt(n)
    // as well as multiple of L1CacheSize
    while (segmentSize < sqrtBound) {
        segmentSize += L1CacheSize;
        numSegments++;
    }
    
    int wheelSize = numSegments * (initialWheelSize + 1);
    std::vector<int> myWheel(wheelSize);
    
    myWheel[0] = 1;
    for (std::size_t i = 0; i < initialWheelSize; i++) {
        myWheel[i+1] = myWheel[i] + wheelDiff30030[i];
    }

    for (std::size_t nSeg = 1; nSeg < numSegments; nSeg++) {
        int wheelIndex = nSeg * (initialWheelSize + 1);
        int wheelOffSet = nSeg * L1CacheSize;
        for (std::size_t i = 0; i <= initialWheelSize; i++)
            myWheel[wheelIndex + i] = myWheel[i] + wheelOffSet;
    }
    
    int miniReserve = (int) floor((double) 2*sqrtBound/log((double)sqrtBound));
    std::vector<typePrime> smallPrimes, nextStrt;
    smallPrimes.reserve(miniReserve);
    nextStrt.reserve(miniReserve);

    uint_fast64_t k, j;
    unsigned long int flrMaxNum = segmentSize * floor((double) n/segmentSize);

    // Create vector of primes up to the sqrt(n)
    // discarding 2, and tacking on the next
    // prime after the sqrt(n)
    lastP = 2;
    if (n < INT_MAX) {
        for (std::size_t i = 0; lastP <= sqrtBound; i++) {
            lastP += primesDiffSeg[i];
            smallPrimes.push_back(lastP);
        }
    } else {
        int sqrtIntMax = floor(sqrt((double) INT_MAX));

        for (std::size_t i = 0; lastP < sqrtIntMax; i++) {
            lastP += primesDiffSeg[i];
            smallPrimes.push_back(lastP);
        }

        // The number, 225, comes from the observation that
        // the largest prime gap less than 100 million is
        // 219 @ 47,326,693. This is important because we
        // need to guarantee that we obtain the smallest
        // prime greater than sqrt(n). We know the that
        // largest number that this algorithm can except
        // is 2^53 - 1, which gives a sqrt of 94,906,266
        int myLimit = sqrtBound + 225;
        std::vector<char> sqrtIsPrime(myLimit + 1, 1);
        lastP = smallPrimes[0];

        for (std::size_t p = 0; lastP * lastP <= myLimit; p++) {
            lastP = smallPrimes[p];
            for (j = (int) (lastP * lastP); j <= myLimit; j += lastP)
                sqrtIsPrime[j] = 0;
        }

        int currentSize = smallPrimes.size();
        lastP = smallPrimes[currentSize - 1];

        for (int j = lastP + 2; j <= myLimit; j += 2)
            if (sqrtIsPrime[j])
                smallPrimes.push_back(j);
    }

    // vector used for sieving
    std::vector<char> sieve(segmentSize, 1);
    sieve[1] = 0;
    upperBnd = segmentSize;

    std::size_t p = 0;
    sqrPrime = smallPrimes[p] * smallPrimes[p];
    typePrime lowerBnd = 0;

    if (minNum > 2) {
        lowerBnd = (typePrime) (segmentSize * floor((double) m/segmentSize));
        upperBnd = std::min(lowerBnd + segmentSize, maxNum);

        typePrime myStart;

        while (sqrPrime <= upperBnd) {
            if (lowerBnd > sqrPrime) {
                double remTest = std::fmod(lowerBnd, smallPrimes[p]);
                if (remTest == 0.0) {
                    myStart = 0;
                } else {
                    myStart = (typePrime) (smallPrimes[p] - remTest);
                }
                if (fmod(lowerBnd + myStart, 2) == 0.0)
                    myStart += smallPrimes[p];
            } else {
                myStart = (typePrime) (sqrPrime - lowerBnd);
            }

            nextStrt.push_back(myStart);
            p++;
            sqrPrime = (typePrime) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = (typePrime) (j - segmentSize);
        }

        std::size_t q = 0;
        while (((typePrime) (myWheel[q]) + lowerBnd) < minNum) {q++;}

        if (upperBnd < flrMaxNum) {
            for (; q < wheelSize; q++)
                if (sieve[myWheel[q]])
                    myPrimes.push_back((typePrime) (myWheel[q]) + lowerBnd);
        } else {
            for (; ((typePrime) (myWheel[q]) + lowerBnd) < maxNum; q++)
                if (sieve[myWheel[q]])
                    myPrimes.push_back((typePrime) (myWheel[q]) + lowerBnd);
        }

        std::fill(sieve.begin(), sieve.end(), 1);
        lowerBnd += segmentSize;
    }

    for (; lowerBnd < flrMaxNum; lowerBnd += segmentSize) {
        // current segment = interval [lowerBnd, upperBnd]
        upperBnd = lowerBnd + segmentSize;

        // sieve the current segment && sieving primes <= sqrt(upperBnd)
        while (sqrPrime <= upperBnd) {
            nextStrt.push_back((typePrime) (sqrPrime - lowerBnd));
            p++;
            sqrPrime = (typePrime) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = (typePrime) (j - segmentSize);
        }

        // Check the sieving interval only at indices that
        // are relatively prime to 30030 = 2*3*5*7*11*13
        for (std::size_t q = 0; q < wheelSize; q++)
            if (sieve[myWheel[q]])
                myPrimes.push_back((typePrime) (myWheel[q]) + lowerBnd);

        std::fill(sieve.begin(), sieve.end(), 1);
    }

    // Get remaining primes that are greater than flrMaxNum and less than maxNum
    if (lowerBnd < maxNum) {
        while (sqrPrime <= maxNum) {
            nextStrt.push_back((typePrime) (sqrPrime - lowerBnd));
            p++;
            sqrPrime = (typePrime) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = (typePrime) (j - segmentSize);
        }

        for (std::size_t q = 0; ((typePrime) (myWheel[q]) + lowerBnd) < maxNum; q++)
            if (sieve[myWheel[q]])
                myPrimes.push_back((typePrime) (myWheel[q]) + lowerBnd);
    }

    return wrap(myPrimes);
}

// [[Rcpp::export]]
List PrimeFactorizationListRcpp (SEXP n) {
    int m;
    double mTest;

    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > INT_MAX) {stop("n must be less than 2^31");}
    
    if (mTest <= 0) {
        stop("n must be positive");
    } else if (mTest < 2) {
        List trivialReturn;
        return trivialReturn;
    }
    
    m = (int)ceil(mTest);
    
    std::vector<std::vector<int> > MyPrimeList(m, std::vector<int>());
    std::vector<std::vector<int> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    IntegerVector primes = AllPrimesCpp<IntegerVector>(1, m);
    IntegerVector::iterator p;
    int i, j, limit, myStep, myMalloc;
    myMalloc = (int)floor(log2((double)m));
    double myLogN = log((double)m);
    
    for (it2d = MyPrimeList.begin(); it2d < itEnd; it2d++) {
        it2d -> reserve(myMalloc);
    }

    for (p = primes.begin(); p < primes.end(); p++) {
        limit = (int)trunc(myLogN/log((double)*p));
        for (i = 1; i <= limit; i++) {
            myStep = (int)pow((double)*p,i);
            for (j = myStep; j <= m; j += myStep) {
                MyPrimeList[j-1].push_back(*p);
            }
        }
    }

    return wrap(MyPrimeList);
}

// [[Rcpp::export]]
IntegerVector EulerPhiSieveRcpp (SEXP n) {
    int m;
    double mTest;

    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > INT_MAX) {stop("n must be less than 2^31");}
    
    if (mTest <= 0) {
        stop("n must be positive");
    } else if (mTest <= 1) {
        IntegerVector trivialReturn(1, 1);
        return trivialReturn;
    }
    
    m = (int)ceil(mTest);
    
    IntegerVector starterSequence = Rcpp::seq(1, m);
    NumericVector EulerPhis = as<NumericVector>(starterSequence);
    std::vector<int> EulerInt(m);
    IntegerVector primes = AllPrimesCpp<IntegerVector>(1, m);
    IntegerVector::iterator p;
    int j;

    for (p = primes.begin(); p < primes.end(); p++) {
        for (j = *p; j <= m; j += *p) {
            EulerPhis[j-1] *= (1.0 - 1.0/(*p));
        }
    }

    for (j = 0; j < m; j++) {EulerInt[j] = round(EulerPhis[j]);}

    return wrap(EulerInt);
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2) {
    double bound1, bound2, myMax, myMin;

    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            bound1 = as<double>(Rb1);
            break;
        }
        case INTSXP: {
            bound1 = as<double>(Rb1);
            break;
        }
        default: {
            stop("bound1 must be of type numeric or integer");
        }
    }

    if (Rf_isNull(Rb2)) {
        if (bound1 <= 0) {stop("n must be positive");}
        
        if (bound1 < 2) {
            IntegerVector v;
            return v;
        } else {
            if (bound1 > INT_MAX)
                return AllPrimesCpp<NumericVector>((double) 1, bound1);
            
            int intBound1 = (int) bound1;
            return AllPrimesCpp<IntegerVector>(1, intBound1);
        }
    } else {
        if (bound1 <= 0 || bound1 > Significand53) {stop("bound1 must be a positive number less than 2^53");}
        
        switch(TYPEOF(Rb2)) {
            case REALSXP: {
                bound2 = as<double>(Rb2);
                break;
            }
            case INTSXP: {
                bound2 = as<double>(Rb2);
                break;
            }
            default: {
                stop("bound2 must be of type numeric or integer");
            }
        }
        if (bound2 <= 0 || bound2 > Significand53) {stop("bound2 must be a positive number less than 2^53");}

        if (bound1 > bound2) {
            myMax = bound1;
            myMin = bound2;
        } else {
            myMax = bound2;
            myMin = bound1;
        }
        
        myMin = ceil(myMin);
        myMax = floor(myMax);

        if (myMax <= 1) {
            IntegerVector z;
            return z;
        }
        
        if (myMax < INT_MAX) {
            if (myMin <= 2) {
                return AllPrimesCpp<IntegerVector>((int) 1, (int) myMax);
            } else {
                return AllPrimesCpp<IntegerVector>((int) myMin, (int) myMax);
            }
        } else {
            if (myMin <= 2) {
                return AllPrimesCpp<NumericVector>((double) 1, myMax);
            } else {
                return AllPrimesCpp<NumericVector>(myMin, myMax);
            }
        }
    }
}
