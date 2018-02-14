#include <Rcpp.h>
#include <math.h>
#include <stdint.h>
#include "PrimesSegSieve.h"
#include "WheelInt30030.h"

using namespace Rcpp;

// "AllPrimesCpp" implements a segmented version of the Sieve of
// Eratosthenes (original implementation authored by Kim Walisch).
// An overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// The official github repo is:
//                      https://github.com/kimwalisch/primesieve

const double Significand53 = 9007199254740991.0;

// 2*3*5*7*11*13...the primes used to generate the wheel
const int L1CacheSize = 30030;

// Number of elements in the static array wheelDiff30030
#define initialWheelSize (sizeof(wheelDiff30030) / sizeof(wheelDiff30030[0]))

template <typename typeRcpp, typename typePrime>
typeRcpp AllPrimesCpp (typePrime m, typePrime n) {

    typePrime sqrtBound = floor(sqrt((double)n));
    typePrime segmentSize = L1CacheSize;
    unsigned int numSegments = 1;
    std::vector<typePrime> myPrimes;
    typePrime minNum = m, maxNum = n;
    typePrime myRange = maxNum - minNum + 1;
    
    // // Percentages obtained here: https://en.wikipedia.org/wiki/Prime-counting_function
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

    int_fast64_t sqrPrime, lastP;
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
    if (segmentSize < sqrtBound) {
        numSegments = (unsigned int) ceil((double) sqrtBound/ (double) segmentSize);
        segmentSize = numSegments * L1CacheSize;
    }

    unsigned int wheelSize = numSegments * (initialWheelSize + 1);
    std::vector<unsigned int> myWheel(wheelSize, 1u);

    for (std::size_t i = 0; i < initialWheelSize; i++)
        myWheel[i+1] = myWheel[i] + wheelDiff30030[i];

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

    int_fast64_t k, j;
    typePrime flrMaxNum = segmentSize * floor((double) n/segmentSize);

    // Create vector of primes up to the sqrt(n)
    // discarding 2, and tacking on the next
    // prime after the sqrt(n)
    lastP = 2;
    if (maxNum < INT_MAX) {
        for (std::size_t i = 0; lastP <= sqrtBound; i++) {
            lastP += primesDiffSeg[i];
            smallPrimes.push_back(lastP);
        }
    } else {
        int_fast64_t sqrtIntMax = floor(sqrt((double) INT_MAX));

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
        int_fast64_t myLimit = sqrtBound + 225;
        std::vector<char> sqrtIsPrime(myLimit + 1, 1);
        lastP = smallPrimes[0];

        for (std::size_t p = 0; lastP * lastP <= myLimit; p++) {
            lastP = smallPrimes[p];
            for (j = (int) (lastP * lastP); j <= myLimit; j += lastP)
                sqrtIsPrime[j] = 0;
        }

        int currentSize = smallPrimes.size();
        lastP = smallPrimes[currentSize - 1];

        for (j = lastP + 2; j <= myLimit; j += 2)
            if (sqrtIsPrime[j])
                smallPrimes.push_back(j);
    }

    // vector used for sieving
    std::vector<char> sieve(segmentSize, 1);
    sieve[1] = 0;

    std::size_t p = 0;
    sqrPrime = (int_fast64_t) (smallPrimes[p] * smallPrimes[p]);
    typePrime lowerBnd = 0, upperBnd = segmentSize;

    if (minNum > 2) {
        lowerBnd = (typePrime) (segmentSize * floor((double) m/segmentSize));
        upperBnd = std::min(lowerBnd + segmentSize, maxNum);

        typePrime myStart;

        while (sqrPrime <= upperBnd) {
            if (lowerBnd > sqrPrime) {
                typePrime remTest = lowerBnd % smallPrimes[p];
                if (remTest == 0) {
                    myStart = 0;
                } else {
                    myStart = (typePrime) (smallPrimes[p] - remTest);
                }
                if (((lowerBnd + myStart) % 2) == 0)
                    myStart += smallPrimes[p];
            } else {
                myStart = (typePrime) (sqrPrime - lowerBnd);
            }

            nextStrt.push_back(myStart);
            p++;
            sqrPrime = (int_fast64_t) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = (int_fast64_t) (smallPrimes[i] * 2); j < segmentSize; j += k)
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
            for (; ((typePrime) (myWheel[q]) + lowerBnd) <= maxNum; q++)
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
            sqrPrime = (int_fast64_t) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = (int_fast64_t) (smallPrimes[i] * 2); j < segmentSize; j += k)
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
            sqrPrime = (int_fast64_t) (smallPrimes[p] * smallPrimes[p]);
        }

        for (std::size_t i = 5; i < p; i++) {
            j = nextStrt[i];
            for (k = (int_fast64_t) (smallPrimes[i] * 2); j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = (typePrime) (j - segmentSize);
        }

        for (std::size_t q = 0; ((typePrime) (myWheel[q]) + lowerBnd) <= maxNum; q++)
            if (sieve[myWheel[q]])
                myPrimes.push_back((typePrime) (myWheel[q]) + lowerBnd);
    }
    
    return wrap(myPrimes);
}

template <typename typeInt>
inline typeInt getStartingIndex (typeInt lowerB,
                          typeInt step, typeInt myPrime) {

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

template <typename typeInt>
List PrimeFactorizationList (typeInt m, typeInt n, bool keepNames) {
    
    typeInt myRange = n;
    myRange += (1 - m);
    
    std::vector<std::vector<typeInt> > MyPrimeList(myRange, std::vector<typeInt>());
    typename std::vector<std::vector<typeInt> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    
    int_fast32_t sqrtBound = floor(sqrt((double)n));
    IntegerVector::iterator p;
    
    typeInt j, limit, myStep, myMalloc, priTypeInt, myNum = m;
    double myLogN = log((double)n);
    unsigned int mySize;
    
    if (n > 1000)
        myMalloc = (typeInt) round(myLogN / log(6.0));
    else
        myMalloc = 4;
    
    std::vector<typeInt> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        for (std::size_t i = 0, j = myNum; j <= n; j++, i++)
            myNames[i] = j;
    }
    
    if (myNum < 2){
        myNum++;
        it2d = MyPrimeList.begin() + 1;
    } else {
        it2d = MyPrimeList.begin();
    }
    
    for (; it2d < itEnd; it2d++, myNum++) {
        it2d -> reserve(myMalloc);
        it2d -> push_back(myNum);
    }
    
    if (n > 3) {
        IntegerVector primes = AllPrimesCpp<IntegerVector>((int_fast32_t) 1, sqrtBound);
    
        if (m < 2) {
            for (p = primes.begin(); p < primes.end(); p++) {
                limit = (typeInt) trunc(myLogN/log((double)*p));
                for (std::size_t i = 1; i <= limit; i++) {
                    myStep = (typeInt) pow((double)*p, i);
                    priTypeInt = (typeInt) (*p);
                    for (j = (myStep - 1); j < n; j += myStep) {
                        mySize = MyPrimeList[j].size() - 1;
                        if (MyPrimeList[j][mySize] > priTypeInt) {
                            MyPrimeList[j][mySize] /= priTypeInt;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, priTypeInt);
                        }
                    }
                }
            }
        } else {
            typeInt myStart;
    
            for (p = primes.begin(); p < primes.end(); p++) {
                limit = (typeInt) trunc(myLogN/log((double)*p));
                
                for (std::size_t i = 1; i <= limit; i++) {
                    myStep = (typeInt) pow((double)*p, i);
                    priTypeInt = (typeInt) (*p);
                    myStart = getStartingIndex(m, myStep, priTypeInt);
    
                    for (j = myStart; j < myRange; j += myStep) {
                        mySize = MyPrimeList[j].size() - 1;
                        if (MyPrimeList[j][mySize] > priTypeInt) {
                            MyPrimeList[j][mySize] /= priTypeInt;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, priTypeInt);
                        }
                    }
                }
            }
        }
    }
    
    Rcpp::List myList = wrap(MyPrimeList);
    if (keepNames)
        myList.attr("names") = myNames;
    
    return myList;
}

template <typename typeRcpp, typename typeInt>
typeRcpp EulerPhiSieveCpp (typeInt m, typeInt n, bool keepNames) {
    
    typeInt myRange = n;
    myRange += (1 - m);
    
    IntegerVector primes;
    IntegerVector::iterator p;
    
    typeInt j, limit, myStep, priTypeInt, myNum = m;
    double myLogN = log((double)n);
    
    std::vector<typeInt> numSeq, EulerPhis(myRange);
    std::vector<typeInt> myNames;
    
    if (keepNames) {
        myNames.resize(myRange);
        for (std::size_t i = 0, j = myNum; j <= n; j++, i++)
            myNames[i] = j;
    }
    
    for (std::size_t i = 0, j = myNum; j <= n; j++, i++)
        EulerPhis[i] = j;
    
    if (m < 2) {
        primes  = AllPrimesCpp<IntegerVector>((typeInt) 1, n);
        
        for (p = primes.begin(); p < primes.end(); p++) {
            myStep = (typeInt) (*p);
            for (j = (myStep - 1); j < n; j += myStep) {
                myNum = EulerPhis[j] / myStep;
                EulerPhis[j] -= myNum;
            }
        }
    } else if (n > 3) {
        numSeq = EulerPhis;
        int_fast32_t sqrtBound = floor(sqrt((double)n));
        primes = AllPrimesCpp<IntegerVector>((int_fast32_t) 1, sqrtBound);
        typeInt myStart;
    
        for (p = primes.begin(); p < primes.end(); p++) {
            limit = (typeInt) trunc(myLogN/log((double)*p));
            myStep = (typeInt) (*p);
            priTypeInt = (typeInt) (*p);
            myStart = getStartingIndex(m, myStep, priTypeInt);
    
            for (j = myStart; j < myRange; j += myStep) {
                numSeq[j] /= priTypeInt;
                myNum = EulerPhis[j] / priTypeInt;
                EulerPhis[j] -= myNum;
            }
    
            for (std::size_t i = 2; i <= limit; i++) {
                myStep = (typeInt) pow((double)*p, i);
                myStart = getStartingIndex(m, myStep, priTypeInt);
    
                for (j = myStart; j < myRange; j += myStep)
                    numSeq[j] /= priTypeInt;
            }
        }

        for (std::size_t i = 0; i < myRange; i++) {
            if (numSeq[i] > 1) {
                myNum = EulerPhis[i] / numSeq[i];
                EulerPhis[i] -= myNum;
            }
        }
    } else { // edge case where m,n = 2 or 3
        for (std::size_t i = 0; i < myRange; i++)
            EulerPhis[i]--;
    }
    
    typeRcpp myVector = wrap(EulerPhis);
    if (keepNames)
        myVector.attr("names") = myNames;
        
    return myVector;
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2, 
                       SEXP RIsList, SEXP RIsEuler, SEXP RNamed) {
    double bound1, bound2, myMax, myMin;
    bool isList = false, isEuler = false, isNamed = false;

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
    
    isList = as<bool>(RIsList);
    isEuler = as<bool>(RIsEuler);
    isNamed = as<bool>(RNamed);

    if (Rf_isNull(Rb2)) {
        if (bound1 <= 0) {stop("n must be positive");}
        myMax = floor(bound1);
        
        if (isList) {
            if (myMax < 2) {
                IntegerVector v;
                return v;
            } else {
                if (bound1 > (INT_MAX - 1)) {
                    return PrimeFactorizationList((int_fast64_t) 1,
                                                  (int_fast64_t) myMax, isNamed);
                }
                return PrimeFactorizationList((int_fast32_t) 1,
                                              (int_fast32_t) myMax, isNamed);
            }
        } else if (isEuler) {
            if (myMax <= 1) {
                IntegerVector trivialReturn(1, 1);
                return trivialReturn;
            } else {
                if (bound1 > (INT_MAX - 1)) {
                    return EulerPhiSieveCpp<NumericVector>((int_fast64_t) 1,
                                                       (int_fast64_t) myMax, isNamed);
                }
                return EulerPhiSieveCpp<IntegerVector>((int_fast32_t) 1,
                                                   (int_fast32_t) myMax, isNamed);
            }
        } else {
            if (myMax < 2) {
                IntegerVector v;
                return v;
            } else {
                if (bound1 > (INT_MAX - 1)) {
                    return AllPrimesCpp<NumericVector>((int_fast64_t) 1,
                                                       (int_fast64_t) myMax);
                }
                return AllPrimesCpp<IntegerVector>((int_fast32_t) 1,
                                                   (int_fast32_t) myMax);
            }
        }
    } else {
        if (bound1 <= 0 || bound1 > Significand53)
            stop("bound1 must be a positive number less than 2^53");
        
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
        if (bound2 <= 0 || bound2 > Significand53)
            stop("bound2 must be a positive number less than 2^53");

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
                std::vector<std::vector<int> > trivialRet(1, std::vector<int>(1, 1));
                Rcpp::List z = wrap(trivialRet);
                if (isNamed)
                    z.attr("names") = 1;
                
                return z;
            } else {
                if (bound1 > (INT_MAX - 1)) {
                    return PrimeFactorizationList((int_fast64_t) myMin,
                                                  (int_fast64_t) myMax, isNamed);
                }
                return PrimeFactorizationList((int_fast32_t) myMin,
                                                   (int_fast32_t) myMax, isNamed);
            }
        } else if (isEuler) {
            if (myMax <= 1) {
                IntegerVector z(1, 1);
                if (isNamed)
                    z.attr("names") = 1;
                
                return z;
            }
            
            if (bound1 > (INT_MAX - 1)) {
                return EulerPhiSieveCpp<NumericVector>((int_fast64_t) myMin,
                                                       (int_fast64_t) myMax, isNamed);
            }
            return EulerPhiSieveCpp<IntegerVector>((int_fast32_t) myMin,
                                                   (int_fast32_t) myMax, isNamed);
        } else {
            if (myMax <= 1) {
                IntegerVector z;
                return z;
            }
            if (myMin == myMax) {myMax++;}
            if (myMin <= 2) myMin = 1;
            if (myMax < (INT_MAX - 1)) {
                    return AllPrimesCpp<IntegerVector>((int_fast32_t) myMin,
                                                       (int_fast32_t) myMax);
            } else {
                    return AllPrimesCpp<NumericVector>((int_fast64_t) myMin, 
                                                       (int_fast64_t) myMax);
            }
        }
    }
}
