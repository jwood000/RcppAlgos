#include <PrimesUtils.h>
#include <GetFacsUtils.h>
#include <CleanConvert.h>
#include <Wheel30030.h>
#include <thread>

// "PrimeSieve" implements a simple segmented version of the Sieve of 
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve

template <typename typePrime>
void PrimeSieveSmall(int_fast64_t minNum, int_fast64_t maxNum,
                     const std::vector<int_fast64_t> &sievePrimes,
                     std::vector<typePrime> &myPrimes, std::vector<int_fast64_t> &nextStrt,
                     int_fast64_t &sqrPrime, unsigned long int &p) {
    
    int_fast64_t segSize = static_cast<int_fast64_t>(L1CacheSize);
    std::size_t numSegs = nWheelsPerSeg, szWheel210 = wheelSize;
    std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum),
                                            static_cast<double>(maxNum));
    myPrimes.reserve(myReserve);
    
    if (maxNum <= smallPrimeBase[lastSmlPri]) {
        std::size_t ind = 0;
        for (; smallPrimeBase[ind] < minNum; ++ind) {}
        
        for (; smallPrimeBase[ind] <= maxNum; ++ind) 
            myPrimes.push_back(static_cast<typePrime>(smallPrimeBase[ind]));
    } else {
        if (minNum < 13) {
            std::size_t ind = 0;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}
            
            for (; smallPrimeBase[ind] < 10; ++ind)
                myPrimes.push_back(static_cast<typePrime>(smallPrimeBase[ind]));
        }
        
        int_fast64_t flrMaxNum = static_cast<int_fast64_t>(segSize * std::floor(static_cast<double>(maxNum) / segSize));
        
        // vector used for sieving
        std::vector<char> sieve(segSize, 1);
        if (minNum < 2) sieve[1] = 0;
        
        int_fast64_t lowerBnd = static_cast<int_fast64_t>(segSize * std::floor(static_cast<double>(minNum) / segSize));
        int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        int_fast64_t myNum = 1 + lowerBnd;

        if (nextStrt.empty() && minNum > 2) {
            int_fast64_t myIndex;
            
            for (; sqrPrime <= upperBnd; ++p) {
                if (lowerBnd > sqrPrime) {
                    int_fast64_t remTest = lowerBnd % sievePrimes[p - 1];
                    if (remTest == 0) {
                        myIndex = sievePrimes[p - 1];
                    } else {
                        myIndex = sievePrimes[p - 1] - remTest;
                        if ((myIndex % 2) == 0) {myIndex += sievePrimes[p - 1];}
                    }
                } else {
                    myIndex = sqrPrime - lowerBnd;
                }
                
                nextStrt.push_back(myIndex);
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            if (upperBnd < flrMaxNum) {
                for (std::size_t q = 0; q < numSegs; ++q)
                    for (std::size_t w = 0; w < szWheel210; myNum += wheel210[w], ++w)
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back(static_cast<typePrime>(myNum));
            } else {
                for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q)
                    for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += wheel210[w], ++w)
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back(static_cast<typePrime>(myNum));
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
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < numSegs; ++q)
                for (std::size_t w = 0; w < szWheel210; myNum += wheel210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back(static_cast<typePrime>(myNum));
            
            std::fill(sieve.begin(), sieve.end(), 1);
        }
        
        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += wheel210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back(static_cast<typePrime>(myNum));
        }
    }
}

template <typename typePrime>
void PrimeSieveMedium(int_fast64_t minNum, int_fast64_t maxNum, 
                      const std::vector<int_fast64_t> &sievePrimes,
                      std::vector<typePrime> &myPrimes, std::vector<int_fast64_t> &nextStrt,
                      int_fast64_t &sqrPrime, unsigned long int &p, std::size_t nSegOne) {
    
    double dblMax = static_cast<double>(maxNum);
    std::size_t nWheels = nSegOne * nWheelsPerSeg, szWheel210 = wheelSize;;
    
    int_fast64_t segSize = static_cast<int_fast64_t>(nSegOne * L1CacheSize);
    std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum), dblMax);
    myPrimes.reserve(myReserve);
    
    int_fast64_t flrMaxNum = segSize * static_cast<int_fast64_t>(std::floor(dblMax / segSize));
    
    // vector used for sieving
    std::vector<bool> sieve(segSize, true);
    
    if (minNum < 13) {
        sieve[1] = false;
        std::size_t ind = 0;
        for (; smallPrimeBase[ind] < minNum; ++ind) {}
        
        for (; smallPrimeBase[ind] < 10; ++ind)
            myPrimes.push_back(static_cast<typePrime>(smallPrimeBase[ind]));
    }
    
    int_fast64_t lowerBnd = static_cast<int_fast64_t>(segSize * std::floor(static_cast<double>(minNum) / segSize));
    int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
    int_fast64_t myNum = 1 + lowerBnd;
    
    if (nextStrt.empty()) {
        int_fast64_t myIndex;
        
        for (; sqrPrime <= upperBnd; ++p) {
            if (lowerBnd > sqrPrime) {
                int_fast64_t remTest = lowerBnd % sievePrimes[p - 1];
                if (remTest == 0) {
                    myIndex = sievePrimes[p - 1];
                } else {
                    myIndex = sievePrimes[p - 1] - remTest;
                    if ((myIndex % 2) == 0) {myIndex += sievePrimes[p - 1];}
                }
            } else {
                myIndex = sqrPrime - lowerBnd;
            }
            
            nextStrt.push_back(myIndex);
            sqrPrime = sievePrimes[p] * sievePrimes[p];
        }
        
        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                sieve[j] = false;
            
            nextStrt[i] = j - segSize;
        }
        
        if (upperBnd < flrMaxNum) {
            for (std::size_t q = 0; q < nWheels; ++q)
                for (std::size_t w = 0; w < szWheel210; myNum += wheel210[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
        } else {
            for (std::size_t q = 0; q < nWheels && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += wheel210[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
        }
        
        std::fill(sieve.begin(), sieve.end(), true);
        lowerBnd += segSize;
    }
    
    for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
        upperBnd = lowerBnd + segSize;
        
        for (; sqrPrime <= upperBnd; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = sievePrimes[p] * sievePrimes[p];
        }
        
        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                sieve[j] = false;
            
            nextStrt[i] = j - segSize;
        }
        
        for (std::size_t q = 0; q < nWheels; ++q)
            for (std::size_t w = 0; w < szWheel210; myNum += wheel210[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back(static_cast<typePrime>(myNum));
        
        std::fill(sieve.begin(), sieve.end(), true);
    }
    
    // Get remaining primes that are greater than flrMaxNum and less than maxNum
    if (lowerBnd < maxNum) {
        for (; sqrPrime <= maxNum; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = sievePrimes[p] * sievePrimes[p];
        }
        
        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                sieve[j] = false;
            
            nextStrt[i] = j - segSize;
        }
        
        for (std::size_t q = 0; q < nWheels && myNum <= maxNum; ++q)
            for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += wheel210[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back(static_cast<typePrime>(myNum));
    }
}

template <typename typePrime>
void PrimeSieveBig(int_fast64_t minNum, int_fast64_t maxNum, 
                   const std::vector<int_fast64_t> &svPriOne,
                   const std::vector<int_fast64_t> &svPriTwo, 
                   std::vector<typePrime> &myPrimes, const int nBigSegs,
                   const std::vector<char> &check30030) {
    
    int_fast64_t segSize = static_cast<int_fast64_t>(nBigSegs * num30030);
    int_fast64_t sz30030 = num30030;
    std::size_t numWheelSegs = nBigSegs, szWheel = szWheel30030;
    std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum),
                                            static_cast<double>(maxNum));
    myPrimes.reserve(myReserve);
    
    int_fast64_t remTest, divTest, myIndex, lowerBnd;
    lowerBnd = static_cast<int_fast64_t>(segSize * (minNum / segSize));
    
    std::size_t svPriOneSize = svPriOne.size();
    std::vector<int_fast64_t> nextStrtOne(svPriOneSize);
    
    for (std::size_t i = 0; i < svPriOneSize; ++i) {
        remTest = lowerBnd % svPriOne[i];
        
        if (remTest == 0) {
            myIndex = svPriOne[i];
        } else {
            myIndex = svPriOne[i] - remTest;
            if (myIndex % 2 == 0) {myIndex += svPriOne[i];}
        }
        
        nextStrtOne[i] = myIndex;
    }
    
    int_fast64_t myRange = (maxNum - lowerBnd) + 1;
    int_fast64_t numCacheSegs = myRange / segSize;
    
    double wholeTest = static_cast<double>(myRange) / segSize;
    
    if (wholeTest != static_cast<int>(wholeTest))
        ++numCacheSegs;
    
    int_fast64_t remPrime, timesTwo;
    int_fast64_t tempInd, maxIndex = myRange + 1;
    bool bKeepGoing;
    
    // Keeps track of which primes will be used in each interval
    std::vector<std::vector<std::size_t>> myBuckets(numCacheSegs,
                                                    std::vector<std::size_t>());

    for (std::size_t i = 0; i < svPriTwo.size(); ++i) {
        remTest = (lowerBnd % svPriTwo[i]);

        if (remTest == 0) {
            myIndex = svPriTwo[i];
        } else {
            myIndex = (svPriTwo[i] - remTest);
            if ((myIndex % 2) == 0) {myIndex += svPriTwo[i];}
        }

        remTest = (myIndex % sz30030) - 1;
        timesTwo = (2 * svPriTwo[i]);
        remPrime = (timesTwo % sz30030);
        bKeepGoing = (myIndex <= maxIndex);
        
        // Populate rest of the buckets
        while (bKeepGoing) {
            while (check30030[remTest] && bKeepGoing) {
                myIndex += timesTwo;
                bKeepGoing = (myIndex <= maxIndex);
                remTest += remPrime;
                if (remTest >= sz30030) {remTest -= sz30030;}
            }
            
            tempInd = myIndex;
            
            if (bKeepGoing) {
                divTest = (myIndex / segSize);
                myIndex -= (divTest * segSize);
                myBuckets[divTest].push_back(static_cast<std::size_t>(myIndex));
            }
            
            myIndex = tempInd + timesTwo;
            bKeepGoing = (myIndex <= maxIndex);
            remTest += remPrime;
            if (remTest >= sz30030) remTest -= sz30030;
        }
    }

    int_fast64_t flrMaxNum = segSize * std::floor(static_cast<double>(maxNum) / segSize);
    int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
    int_fast64_t myNum = 1 + lowerBnd;

    // vector used for sieving
    std::vector<bool> sieve(segSize, true);
    std::size_t strt = 0;

    if (minNum != lowerBnd) {
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;

            nextStrtOne[i] = (j - segSize);
        }

        for (const auto &elem: myBuckets[0])
            sieve[elem] = false;

        if (upperBnd < flrMaxNum) {
            for (std::size_t q = 0; q < numWheelSegs; ++q)
                for (std::size_t w = 0; w < szWheel; myNum += wheel30030[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
        } else {
            for (std::size_t q = 0; q < numWheelSegs && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel && myNum <= maxNum; myNum += wheel30030[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
        }

        std::fill(sieve.begin(), sieve.end(), true);
        lowerBnd += segSize;
        ++strt;
    }

    for (std::size_t v = strt; lowerBnd < flrMaxNum; ++v, lowerBnd += segSize) {
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;

            nextStrtOne[i] = (j - segSize);
        }

        for (const auto &elem: myBuckets[v])
            sieve[elem] = false;
        
        for (std::size_t q = 0; q < numWheelSegs; ++q)
            for (std::size_t w = 0; w < szWheel; myNum += wheel30030[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back(static_cast<typePrime>(myNum));

        std::fill(sieve.begin(), sieve.end(), true);
    }
    
    // Get remaining primes that are greater than flrMaxNum and less than maxNum
    if (lowerBnd < maxNum) {
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;
        }

        for (const auto &elem: myBuckets[numCacheSegs - 1])
            sieve[elem] = false;

        for (std::size_t q = 0; (myNum <= maxNum) && (q < numWheelSegs); ++q)
            for (std::size_t w = 0; (myNum <= maxNum) && (w < szWheel); myNum += wheel30030[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back(static_cast<typePrime>(myNum));
    }
}

void sqrtSmallPrimes(int sqrtBound, std::vector<int_fast64_t> &sievePrimes) {
    std::size_t ind = 1;
    for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
        sievePrimes.push_back(smallPrimeBase[ind]);

    sievePrimes.push_back(smallPrimeBase[ind]);
}

void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                   bool bAddTwo, std::vector<int_fast64_t> &sievePrimes) {

    if (sqrtBound < smallPrimeBase[lastSmlPri]) {
        if (bAddZero) sievePrimes.push_back(0);
        unsigned long int ind = (bAddTwo) ? 0 : 1;
        
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            sievePrimes.push_back(smallPrimeBase[ind]);
        
        if (bAddExtraPrime)
            sievePrimes.push_back(smallPrimeBase[ind]);
    } else {
        int sqrtSqrtBound = static_cast<int>(std::sqrt(sqrtBound));
        std::vector<int_fast64_t> sqrtSievePrimes, sqrtNextStrt;
        sqrtSmallPrimes(sqrtSqrtBound, sqrtSievePrimes);
        
        // The number, 225, comes from the observation that the largest prime
        // gap less than 100 million is 219 @ 47,326,693. This is important
        // because we need to guarantee that we obtain the smallest prime
        // greater than sqrt(n). We know that the largest number that this
        // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266
        
        int_fast64_t myLower = 3, myUpper = sqrtBound;
        if (bAddExtraPrime) {myUpper += 225;}
        std::size_t sqrtReserve = EstimatePiPrime(1.0, static_cast<double>(myUpper));
        sievePrimes.reserve(sqrtReserve);
        
        if (bAddZero) {sievePrimes.push_back(0);}
        if (bAddTwo) {myLower = 1;}
        int_fast64_t sqrtSqrPrime = 9;
        unsigned long int p = 1;
        PrimeSieveSmall(myLower, myUpper, sqrtSievePrimes, 
                        sievePrimes, sqrtNextStrt, sqrtSqrPrime, p);
    }
}

template <typename typePrime>
void MasterPrime(int_fast64_t minNum, int_fast64_t maxNum, int_fast64_t range,
                 int_fast64_t lowerBnd, int_fast64_t upperBnd, int_fast64_t chunkSize, 
                 const std::vector<int_fast64_t> &svPriMain, const std::vector<int_fast64_t> &svPriOne,
                 const std::vector<int_fast64_t> &svPriTwo, const std::vector<char> &check30030,
                 std::vector<typePrime> &primes, std::vector<std::vector<typePrime>> &primeList,
                 int numSects, int nCacheL1, int_fast64_t medCut, int nThreads = 0) {
    
    unsigned long int p = 1, ind = 0;
    int_fast64_t sqrPrime = 9;
    std::vector<int_fast64_t> nextStrt;
    
    if ((range > maxPriPer) && (maxNum > medCut || maxNum <= smallCut)) {
        for (; ind < (numSects - 1); lowerBnd = upperBnd, upperBnd += chunkSize, ++ind) {
            if (upperBnd > medCut)
                PrimeSieveBig(lowerBnd, upperBnd, svPriOne, svPriTwo, primeList[ind], nCacheL1, check30030);
            else if (upperBnd > smallCut)
                PrimeSieveMedium(lowerBnd, upperBnd, svPriMain, primeList[ind], nextStrt, sqrPrime, p, nCacheL1);
            else
                PrimeSieveSmall(lowerBnd, upperBnd, svPriMain, primeList[ind], nextStrt, sqrPrime, p);
        }

        if (maxNum > medCut)
            PrimeSieveBig(lowerBnd, maxNum, svPriOne, svPriTwo, primeList[ind], nCacheL1, check30030);
        else if (upperBnd > smallCut)
            PrimeSieveMedium(lowerBnd, maxNum, svPriMain, primeList[ind], nextStrt, sqrPrime, p, nCacheL1);
        else
            PrimeSieveSmall(lowerBnd, maxNum, svPriMain, primeList[ind], nextStrt, sqrPrime, p);

    } else {
        if (maxNum > medCut)
            PrimeSieveBig(minNum, maxNum, svPriOne, svPriTwo, primes, nCacheL1, check30030);
        else if (maxNum > smallCut)
            PrimeSieveMedium(minNum, maxNum, svPriMain, primes, nextStrt, sqrPrime, p, nCacheL1);
        else
            PrimeSieveSmall(minNum, maxNum, svPriMain, primes, nextStrt, sqrPrime, p);
    }
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp2 (SEXP Rb1, SEXP Rb2, SEXP RNumThreads) {
    double bound1, bound2, myMax, myMin;
    bool Parallel = false;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1 must be of type numeric or integer", false);
    
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
    
    myMin = std::ceil(myMin);
    myMax = std::floor(myMax);
        
    if (myMax <= 1)
        return Rcpp::IntegerVector();
    
    if (myMin <= 2) myMin = 1;
    if (myMin == myMax) {++myMax;}
    
    int numThreads;
    int_fast64_t myRange = static_cast<int_fast64_t>(myMax - myMin);
    int sqrtMedCut = 6 * 4096 * 1024;
    
    int totalThreads = static_cast<int>(std::thread::hardware_concurrency());
    
    int sqrtBound = std::sqrt(myMax);
    std::size_t nCacheL1 = 0u;
    
    if ((myRange < 200000) || (totalThreads < 2)) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        Parallel = true;
        CleanConvert::convertPrimitive(RNumThreads, numThreads, "nThreads must be of type numeric or integer");
        if (numThreads > totalThreads) {numThreads = totalThreads;}
        if (numThreads < 2) {Parallel = false;}
        if (Parallel) {sqrtMedCut /= numThreads;}
    }
    
    int_fast64_t medCut = sqrtMedCut * sqrtMedCut;
    std::vector<int_fast64_t> sievePrimes, sievePrimesOne, sievePrimesTwo;
    std::vector<char> check30030(num30030, 1);
    int_fast64_t mySegSize = static_cast<int_fast64_t>(L1CacheSize);
    
    if (myMax > medCut) {
        sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
        nCacheL1 = maxNumL1Cache;
        mySegSize = num30030;
        
        std::size_t ind = 0u;
        std::size_t limitOne = static_cast<std::size_t>(nCacheL1 * num30030);
        
        // Get the primes that are guaranteed to mark an
        // index in the the every segment interval
        for (; (2 * sievePrimes[ind]) < limitOne; ++ind)
            sievePrimesOne.push_back(sievePrimes[ind]);
        
        // Get the rest
        for (; ind < sievePrimes.size(); ++ind)
            sievePrimesTwo.push_back(sievePrimes[ind]);
        
        for (std::size_t i = 0, ind = 0; i < szWheel30030; ind += wheel30030[i], ++i)
            check30030[ind] = 0;
        
    } else if (myMax > smallCut) {
        sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
        nCacheL1 = static_cast<std::size_t>(std::ceil(
            static_cast<double>(sqrtBound) / L1CacheSize));
    } else {
        sqrtSmallPrimes(sqrtBound, sievePrimes);
    }
    
    int_fast64_t int64Min = static_cast<int_fast64_t>(myMin);
    int_fast64_t int64Max = static_cast<int_fast64_t>(myMax);
    int_fast64_t lowerBnd, upperBnd, tempLower;
    int_fast64_t segSize = static_cast<int_fast64_t>(nCacheL1 * mySegSize);
    
    tempLower = static_cast<int_fast64_t>(segSize * std::floor(myMin / segSize));
    int_fast64_t chunkSize = static_cast<int_fast64_t>(segSize * 
                    std::floor(static_cast<double>(maxPriPer) / segSize));
    
    upperBnd = tempLower + chunkSize;
    lowerBnd = int64Min;
    
    int_fast64_t tempSects = (int64Max - upperBnd) / chunkSize;
    unsigned long int numSects = tempSects + 1;
    std::size_t numPrimes = 0;
    std::vector<unsigned long int> runningCount;
    runningCount.push_back(0u);
    
    if (Parallel) {
        
    }
    
    if (myMax > INT_MAX) {
        std::vector<std::vector<double>> primeList(numSects, std::vector<double>());
        std::vector<double> tempPrime;

        MasterPrime(int64Min, int64Max, myRange, lowerBnd, upperBnd, chunkSize,
                    sievePrimes, sievePrimesOne, sievePrimesTwo, check30030,
                    tempPrime, primeList, numSects, nCacheL1, medCut);

        if ((myRange > maxPriPer) && (myMax > medCut || myMax <= smallCut)) {
            for (int i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }

            Rcpp::NumericVector primes(numPrimes);
            Rcpp::NumericVector::iterator priBeg = primes.begin();

            for (int i = 0; i < numSects; ++i)
                std::copy(primeList[i].begin(), primeList[i].end(), priBeg + runningCount[i]);

            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    } else {
        std::vector<std::vector<int_fast32_t>> primeList(numSects, std::vector<int_fast32_t>());
        std::vector<int_fast32_t> tempPrime;

        MasterPrime(int64Min, int64Max, myRange, lowerBnd, upperBnd, chunkSize,
                    sievePrimes, sievePrimesOne, sievePrimesTwo, check30030,
                    tempPrime, primeList, numSects, nCacheL1, medCut);

        if ((myRange > maxPriPer) && (myMax > medCut || myMax <= smallCut)) {
            for (int i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }

            Rcpp::IntegerVector primes(numPrimes);
            Rcpp::IntegerVector::iterator priBeg = primes.begin();

            for (int i = 0; i < numSects; ++i)
                std::copy(primeList[i].begin(), primeList[i].end(), priBeg + runningCount[i]);

            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    }
    // std::vector<int_fast64_t> sievePrimes;
    // 
    // // if (myMax >= smallLimit) {
    // //     sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
    // //     std::vector<int_fast64_t> primes, sievePrimesOne, sievePrimesTwo;
    // //     std::size_t ind = 0, limitOne = numSegs * segmentSize;
    // //     
    // //     // Get the primes that are guaranteed to mark an
    // //     // index in the the every segment interval
    // //     for ( ; (2 * sievePrimes[ind]) < limitOne; ++ind)
    // //         sievePrimesOne.push_back(sievePrimes[ind]);
    // //     
    // //     // Get the rest
    // //     for ( ; ind < sievePrimes.size(); ++ind)
    // //         sievePrimesTwo.push_back(sievePrimes[ind]);
    // //     
    // // 
    // //     if (Parallel) {
    // //         std::vector<std::thread> myThreads;
    // //         std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    // //         int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
    // //         int_fast64_t upperBnd = myMin + stepSize - 1;
    // // 
    // //         for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    // //             myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    // //                                    upperBnd, std::ref(sievePrimes),
    // //                                    sqrtBound, std::ref(primeList[j]));
    // //         }
    // // 
    // //         myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    // //                                static_cast<int_fast64_t>(myMax, std::ref(sievePrimes),
    // //                                sqrtBound, std::ref(primeList[numThreads - 1]));
    // // 
    // //         for (auto& thr: myThreads)
    // //             thr.join();
    // //         
    // //         std::size_t numPrimes = 0, count = 0;
    // //         std::vector<unsigned long int> sectionSize(numThreads);
    // //         
    // //         for (int i = 0; i < numThreads; ++i) {
    // //             numPrimes += primeList[i].size();
    // //             sectionSize[i] = primeList[i].size();
    // //         }
    // // 
    // //         if (numPrimes < INT_MAX) {
    // //             Rcpp::NumericVector primes(numPrimes);
    // // 
    // //             for (int i = 0; i < numThreads; ++i)
    // //                 for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    // //                     primes[count] = primeList[i][j];
    // // 
    // //             return primes;
    // //         } else {
    // //             return Rcpp::wrap(primeList);
    // //         }
    // //     } else {
    // //         PrimeSieveBig(static_cast<int_fast64_t>(myMin, static_cast<int_fast64_t>(myMax, sievePrimesOne, sievePrimesTwo, primes);
    // //         return Rcpp::wrap(primes);
    // //     }
    // // }
    // //     std::vector<int_fast16_t> small16Primes = sqrt16BasePrimes(sqrtBound);
    // //     std::vector<std::vector<int_fast32_t> > primeList(numThreads, std::vector<int_fast32_t>());
    // //     int_fast32_t lowerBnd = myMin, stepSize = myRange / numThreads;
    // //     int_fast32_t upperBnd = myMin + stepSize - 1;
    // // 
    //     // for (int j = 0; j < (numThreads - 1); ++j) {
    //     //     myThreads.emplace_back(PrimeSieveSmall, lowerBnd, upperBnd,
    //     //                            std::ref(small16Primes), sqrtBound, std::ref(primeList[j]));
    //     //     lowerBnd += stepSize;
    //     //     upperBnd += stepSize;
    //     // }
    //     // 
    //     // myThreads.emplace_back(PrimeSieveSmall, lowerBnd, static_cast<int_fast32_t>(myMax,
    //     //                        std::ref(small16Primes), sqrtBound, std::ref(primeList[numThreads - 1]));
    //     // 
    //     // for (auto& thr: myThreads)
    //     //     thr.join();
    //     // 
    //     // for (int i = 0; i < numThreads; ++i) {
    //     //     numPrimes += primeList[i].size();
    //     //     sectionSize[i] = primeList[i].size();
    //     // }
    //     // 
    //     // Rcpp::IntegerVector primes(numPrimes);
    //     // 
    //     // for (int i = 0; i < numThreads; ++i)
    //     //     for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //     //         primes[count] = primeList[i][j];
    //     // 
    //     // return primes;
    // // }
    // 
    // if (bigP) {
    //     Rcpp::print(Rcpp::wrap("brah"));
    //     sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
    //     std::vector<int_fast64_t> sievePrimesOne, sievePrimesTwo;
    //     std::size_t ind = 0, limitOne = (nCacheL1 * segmentSize);
    //     
    //     // Get the primes that are guaranteed to mark an
    //     // index in the the every segment interval
    //     for (; (2 * sievePrimes[ind]) < limitOne; ++ind)
    //         sievePrimesOne.push_back(sievePrimes[ind]);
    //     
    //     // Get the rest
    //     for (; ind < sievePrimes.size(); ++ind)
    //         sievePrimesTwo.push_back(sievePrimes[ind]);
    //     
    //     if (Parallel) {
    //         if (myRange > 201277440) {
    //             myRange = 201277440;
    //         }
    //         // 
    //         int_fast64_t tempMax = myMin + myRange;
    //         std::vector<Rcpp::NumericVector> primeTemp;
    //         int_fast64_t stepSize = myRange / numThreads;
    //         int_fast64_t lowerBnd = myMin;
    //         int_fast64_t upperBnd = lowerBnd + stepSize - 1;
    //         
    //         Rcpp::print(Rcpp::wrap(tempMax));
    //         Rcpp::print(Rcpp::wrap(lowerBnd));
    //         Rcpp::print(Rcpp::wrap(upperBnd));
    //         
    //         Rcpp::print(Rcpp::wrap(myMax));
    //         
    //         for (; tempMax <= static_cast<int_fast64_t>(myMax); tempMax += myRange) {
    //             // Rcpp::print(Rcpp::wrap(tempMax));
    //             // Rcpp::print(Rcpp::wrap(lowerBnd));
    //             std::vector<std::thread> myThreads;
    //             std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    //         
    //             for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //                 myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                        upperBnd, std::ref(sievePrimesOne),
    //                                        std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
    //             }
    // 
    //             myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                    static_cast<int_fast64_t>(tempMax), std::ref(sievePrimesOne),
    //                                    std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);
    // 
    //             for (auto& thr: myThreads)
    //                 thr.join();
    // 
    //             std::size_t numPrimes = 0, count = 0;
    //             std::vector<unsigned long int> sectionSize(numThreads);
    // 
    //             for (int i = 0; i < numThreads; ++i) {
    //                 numPrimes += primeList[i].size();
    //                 sectionSize[i] = primeList[i].size();
    //             }
    //             
    //             Rcpp::NumericVector primes(numPrimes);
    // 
    //             for (int i = 0; i < numThreads; ++i)
    //                 for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //                     primes[count] = primeList[i][j];
    // 
    //             primeTemp.push_back(primes);
    // 
    //             lowerBnd = tempMax;
    //             upperBnd = lowerBnd + stepSize - 1;
    //         }
    // 
    //         Rcpp::print(Rcpp::wrap("lowerBnd"));
    //         Rcpp::print(Rcpp::wrap(lowerBnd));
    //         Rcpp::print(Rcpp::wrap(tempMax));
    // 
    //         stepSize = (myMax - (tempMax - myRange)) / numThreads;
    //         upperBnd = lowerBnd + stepSize - 1;
    //         std::vector<std::thread> myThreads;
    //         std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    // 
    //         Rcpp::print(Rcpp::wrap(upperBnd));
    //         Rcpp::print(Rcpp::wrap(stepSize));
    // 
    //         for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //             myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                    upperBnd, std::ref(sievePrimesOne),
    //                                    std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
    //         }
    // 
    //         myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                static_cast<int_fast64_t>(myMax), std::ref(sievePrimesOne),
    //                                std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);
    // 
    //         for (auto& thr: myThreads)
    //             thr.join();
    // 
    //         std::size_t numPrimes = 0, count = 0;
    //         std::vector<unsigned long int> sectionSize(numThreads);
    // 
    //         for (int i = 0; i < numThreads; ++i) {
    //             numPrimes += primeList[i].size();
    //             sectionSize[i] = primeList[i].size();
    //         }
    // 
    //         Rcpp::NumericVector primes(numPrimes);
    // 
    //         for (int i = 0; i < numThreads; ++i)
    //             for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //                 primes[count] = primeList[i][j];
    // 
    //         primeTemp.push_back(primes);
    //         
    //         return Rcpp::wrap(primeTemp);
    //         
    //         // int_fast64_t tempMax = myMin + myRange;
    //         // std::vector<Rcpp::NumericVector> primeTemp;
    //         // 
    //         // 
    //         // Rcpp::print(Rcpp::wrap("brah2"));
    //         // std::vector<std::thread> myThreads;
    //         // std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    //         // int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
    //         // int_fast64_t upperBnd = lowerBnd + stepSize - 1;
    //         // 
    //         // for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //         //     myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //         //                            upperBnd, std::ref(sievePrimesOne),
    //         //                            std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
    //         // }
    //         // 
    //         // myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //         //                        static_cast<int_fast64_t>(myMax, std::ref(sievePrimesOne),
    //         //                        std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);
    //         // 
    //         // for (auto& thr: myThreads)
    //         //     thr.join();
    //         // 
    //         // std::size_t numPrimes = 0, count = 0;
    //         // std::vector<unsigned long int> sectionSize(numThreads);
    //         // 
    //         // for (int i = 0; i < numThreads; ++i) {
    //         //     numPrimes += primeList[i].size();
    //         //     sectionSize[i] = primeList[i].size();
    //         // }
    //         // 
    //         // if (numPrimes < INT_MAX) {
    //         //     Rcpp::NumericVector primes(numPrimes);
    //         // 
    //         //     for (int i = 0; i < numThreads; ++i)
    //         //         for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //         //             primes[count] = primeList[i][j];
    //         // 
    //         //     return primes;
    //         // } else {
    //         //     return Rcpp::wrap(primeList);
    //         // }
    //     } else {
    //         Rcpp::print(Rcpp::wrap("brah3"));
    //         std::vector<double> primes;
    //         PrimeSieveBig(static_cast<int_fast64_t>(myMin), static_cast<int_fast64_t>(myMax),
    //                       sievePrimesOne, sievePrimesTwo, primes, nCacheL1);
    //         return Rcpp::wrap(primes);
    //     }
    // } else if (medP) {
    //     Rcpp::print(Rcpp::wrap("brah20"));
    //     sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
    //     
    //     if (Parallel) {
    //         Rcpp::print(Rcpp::wrap("brah22"));
    //         
    //         if (myRange > 1000000000) {
    //             myRange = 1000000000;
    //         }
    //         
    //         int_fast64_t tempMax = myMin + myRange;
    //         std::vector<Rcpp::NumericVector> primeTemp;
    //         int_fast64_t stepSize = myRange / numThreads;
    //         int_fast64_t lowerBnd = myMin;
    //         int_fast64_t upperBnd = lowerBnd + stepSize - 1;
    //         
    //         for (; tempMax < static_cast<int_fast64_t>(myMax); tempMax += myRange) {
    //             Rcpp::print(Rcpp::wrap(tempMax));
    //             Rcpp::print(Rcpp::wrap(lowerBnd));
    //             std::vector<std::thread> myThreads;
    //             std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    //             
    //             for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //                 myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
    //                                        upperBnd, std::ref(sievePrimes),
    //                                        sqrtBound, std::ref(primeList[j]), nCacheL1);
    //             }
    //             
    //             myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
    //                                    static_cast<int_fast64_t>(tempMax), std::ref(sievePrimes),
    //                                    sqrtBound, std::ref(primeList[numThreads - 1]), nCacheL1);
    //             
    //             for (auto& thr: myThreads)
    //                 thr.join();
    //             
    //             std::size_t numPrimes = 0, count = 0;
    //             std::vector<unsigned long int> sectionSize(numThreads);
    //             
    //             for (int i = 0; i < numThreads; ++i) {
    //                 numPrimes += primeList[i].size();
    //                 sectionSize[i] = primeList[i].size();
    //             }
    //             
    //             Rcpp::NumericVector primes(numPrimes);
    //             
    //             for (int i = 0; i < numThreads; ++i)
    //                 for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //                     primes[count] = primeList[i][j];
    //             
    //             primeTemp.push_back(primes);
    //             
    //             lowerBnd = tempMax;
    //             upperBnd = lowerBnd + stepSize - 1;
    //         }
    //         
    //         Rcpp::print(Rcpp::wrap("lowerBnd"));
    //         Rcpp::print(Rcpp::wrap(lowerBnd));
    //         Rcpp::print(Rcpp::wrap(tempMax));
    //         
    //         stepSize = (myMax - (tempMax - myRange)) / numThreads;
    //         upperBnd = lowerBnd + stepSize - 1;
    //         std::vector<std::thread> myThreads;
    //         std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    //         
    //         Rcpp::print(Rcpp::wrap(upperBnd));
    //         Rcpp::print(Rcpp::wrap(stepSize));
    // 
    //         for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //             myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
    //                                    upperBnd, std::ref(sievePrimes),
    //                                    sqrtBound, std::ref(primeList[j]), nCacheL1);
    //         }
    // 
    //         myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
    //                                static_cast<int_fast64_t>(myMax), std::ref(sievePrimes),
    //                                sqrtBound, std::ref(primeList[numThreads - 1]), nCacheL1);
    // 
    //         for (auto& thr: myThreads)
    //             thr.join();
    // 
    //         std::size_t numPrimes = 0, count = 0;
    //         std::vector<unsigned long int> sectionSize(numThreads);
    // 
    //         for (int i = 0; i < numThreads; ++i) {
    //             numPrimes += primeList[i].size();
    //             sectionSize[i] = primeList[i].size();
    //         }
    // 
    //         Rcpp::NumericVector primes(numPrimes);
    // 
    //         for (int i = 0; i < numThreads; ++i)
    //             for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //                 primes[count] = primeList[i][j];
    // 
    //         primeTemp.push_back(primes);
    //         
    //         return Rcpp::wrap(primeTemp);
    //         // if (numPrimes < INT_MAX) {
    //             // Rcpp::NumericVector primes(numPrimes);
    //             // 
    //             // for (int i = 0; i < numThreads; ++i)
    //             //     for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //             //         primes[count] = primeList[i][j];
    //             // 
    //             // return primes;
    //         // } else {
    //         //     return Rcpp::wrap(primeList);
    //         // }
    //     } else {
    //         Rcpp::print(Rcpp::wrap("brah23"));
    //         std::vector<double> primes;
    //         PrimeSieveMedium(static_cast<int_fast64_t>(myMin), static_cast<int_fast64_t>(myMax),
    //                       sievePrimes, sqrtBound, primes, nCacheL1);
    //         return Rcpp::wrap(primes);
    //     }
    // }
    // 
    // sqrtSmallPrimes(sqrtBound, sievePrimes);
    // std::vector<int_fast32_t> primes;
    // 
    // if (Parallel) {
    //     std::vector<std::thread> myThreads;
    //     std::vector<std::vector<int_fast32_t> > primeList(numThreads, std::vector<int_fast32_t>());
    //     int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
    //     int_fast64_t upperBnd = myMin + stepSize - 1;
    //     
    //     for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //         myThreads.emplace_back(PrimeSieveSmall, lowerBnd, upperBnd,
    //                                std::ref(sievePrimes), sqrtBound, std::ref(primeList[j]));
    //     }
    //     
    //     myThreads.emplace_back(PrimeSieveSmall, lowerBnd, static_cast<int_fast64_t>(myMax),
    //                            std::ref(sievePrimes), sqrtBound, std::ref(primeList[numThreads - 1]));
    //     
    //     for (auto& thr: myThreads)
    //         thr.join();
    //     
    //     std::size_t numPrimes = 0, count = 0;
    //     std::vector<unsigned long int> sectionSize(numThreads);
    //     
    //     for (int i = 0; i < numThreads; ++i) {
    //         numPrimes += primeList[i].size();
    //         sectionSize[i] = primeList[i].size();
    //     }
    //     
    //     Rcpp::IntegerVector primes(numPrimes);
    //     
    //     for (int i = 0; i < numThreads; ++i)
    //         for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //             primes[count] = primeList[i][j];
    //     
    //     return primes;
    // } else {
    //     PrimeSieveSmall(static_cast<int_fast64_t>(myMin), static_cast<int_fast64_t>(myMax), sievePrimes, sqrtBound, primes);
    //     return Rcpp::wrap(primes);
    // }
}
