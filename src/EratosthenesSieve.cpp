#include <PrimesUtils.h>
#include <GetFacsUtils.h>
#include <CleanConvert.h>
#include <thread>

// "PrimeSieve" implements a simple segmented version of the Sieve of 
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve

void PrimeSieveSmall(int_fast64_t minNum, int_fast64_t maxNum,
                     const std::vector<int_fast64_t> &sievePrimes,
                     int sqrtBound, std::vector<int_fast32_t> &myPrimes) {
    
    int_fast64_t segSize = segmentSize;
    unsigned long int numSegs = nWheelsPerSeg;
    std::size_t myReserve = EstimatePiPrime((double) minNum, (double) maxNum);
    myPrimes.reserve(myReserve);
    
    if (maxNum <= smallPrimeBase[smlPriBsSize - 1]) {
        std::size_t ind = 0;
        for (; smallPrimeBase[ind] < minNum; ++ind) {}
        
        for (; smallPrimeBase[ind] <= maxNum; ++ind) 
            myPrimes.push_back((int_fast32_t) smallPrimeBase[ind]);
    } else {
        if (minNum < 13) {
            std::size_t ind = 0;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}
            
            for (; smallPrimeBase[ind] < 10; ++ind)
                myPrimes.push_back((int_fast32_t) smallPrimeBase[ind]);
        }
        
        // Ensure segSize is greater than sqrt(n)
        // as well as multiple of L1CacheSize
        if (segSize < sqrtBound) {
            numSegs = (unsigned long int) ceil((double) sqrtBound / sz210);
            segSize = numSegs * sz210;
        }
        
        std::vector<int_fast64_t> nextStrt;
        int_fast64_t flrMaxNum = segSize * std::floor((double) maxNum / segSize);
        
        // vector used for sieving
        std::vector<char> sieve(segSize, 1);
        sieve[1] = 0;
        
        int_fast64_t sqrPrime = (int_fast64_t) (3 * 3);
        int_fast64_t lowerBnd = (int_fast64_t) (segSize * std::floor((double) minNum / segSize));
        int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        int_fast64_t myNum = 1 + lowerBnd;
        int p = 1;
        
        if (minNum > 2) {
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
                for (std::size_t q = 0; q < numSegs; ++q) {
                    for (std::size_t w = 0; w < wheelSize; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((int_fast32_t) myNum);
                            
                        myNum += wheel210[w];
                    }
                }
            } else {
                for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                    for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((int_fast32_t) myNum);
                            
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
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                int_fast64_t j = nextStrt[i];
                for (int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < numSegs; ++q) {
                for (std::size_t w = 0; w < wheelSize; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((int_fast32_t) myNum);
                    
                    myNum += wheel210[w];
                }
            }
            
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
            
            for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((int_fast32_t) myNum);
                    
                    myNum += wheel210[w];
                }
            }
        }
    }
}

template <typename typePrime>
void PrimeSieveMedium(int_fast64_t minNum, int_fast64_t maxNum,
                      const std::vector<int_fast64_t> &sievePrimes,
                      int sqrtBound, std::vector<typePrime> &myPrimes, int nCacheL1) {
    
    int_fast64_t segSize = nCacheL1 * segmentSize;
    unsigned long int numSegs = nCacheL1 * nWheelsPerSeg;
    std::size_t myReserve = EstimatePiPrime((double) minNum, (double) maxNum);
    myPrimes.reserve(myReserve);
    
    if (maxNum <= smallPrimeBase[smlPriBsSize - 1]) {
        std::size_t ind = 0;
        for (; smallPrimeBase[ind] < minNum; ++ind) {}
        
        for (; smallPrimeBase[ind] <= maxNum; ++ind) 
            myPrimes.push_back((typePrime) smallPrimeBase[ind]);
    } else {
        if (minNum < 13) {
            std::size_t ind = 0;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}
            
            for (; smallPrimeBase[ind] < 10; ++ind)
                myPrimes.push_back((typePrime) smallPrimeBase[ind]);
        }
        
        // Ensure segSize is greater than sqrt(n)
        // as well as multiple of L1CacheSize
        if (segSize < sqrtBound) {
            numSegs = (unsigned long int) ceil((double) sqrtBound / sz210);
            segSize = numSegs * sz210;
        }
        
        std::vector<int_fast64_t> nextStrt;
        int_fast64_t flrMaxNum = segSize * std::floor((double) maxNum / segSize);
        
        // vector used for sieving
        std::vector<bool> sieve(segSize, true);
        sieve[1] = false;
        
        int_fast64_t sqrPrime = (int_fast64_t) (3 * 3);
        int_fast64_t lowerBnd = (int_fast64_t) (segSize * std::floor((double) minNum / segSize));
        int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        int_fast64_t myNum = 1 + lowerBnd;
        int p = 1;
        
        if (minNum > 2) {
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
                for (std::size_t q = 0; q < numSegs; ++q) {
                    for (std::size_t w = 0; w < wheelSize; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((typePrime) myNum);
                            
                            myNum += wheel210[w];
                    }
                }
            } else {
                for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                    for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                        if (myNum >= minNum)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back((typePrime) myNum);
                            
                            myNum += wheel210[w];
                    }
                }
            }
            
            std::fill(sieve.begin(), sieve.end(), true);
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
                    sieve[j] = false;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < numSegs; ++q) {
                for (std::size_t w = 0; w < wheelSize; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((typePrime) myNum);
                    
                    myNum += wheel210[w];
                }
            }
            
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
            
            for (std::size_t q = 0; q < numSegs && myNum <= maxNum; ++q) {
                for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back((typePrime) myNum);
                    
                    myNum += wheel210[w];
                }
            }
        }
    }
}

void PrimeSieveBig(int_fast64_t minNum, int_fast64_t maxNum, const std::vector<int_fast64_t> &svPriOne,
                   const std::vector<int_fast64_t> &svPriTwo, std::vector<double> &myPrimes, int nCacheL1) {
    
    int_fast64_t segSize = (nCacheL1 * segmentSize);
    unsigned long int numWheelSegs = (nCacheL1 * nWheelsPerSeg);
    std::size_t myReserve = EstimatePiPrime((double) minNum, (double) maxNum);
    myPrimes.reserve(myReserve);
    
    if (minNum < 13) {
        std::size_t ind = 0;
        for (; svPriOne[ind] < minNum; ++ind) {}
        
        for (; svPriOne[ind] < 10; ++ind)
            myPrimes.push_back((double) svPriOne[ind]);
    }
    
    int_fast64_t remTest, divTest, myIndex, lowerBnd;
    lowerBnd = (int_fast64_t) segSize * (minNum / segSize);
    
    unsigned long int svPriOneSize = svPriOne.size();
    std::vector<int_fast64_t> nextStrtOne(svPriOneSize);
    std::vector<int_fast64_t> nextStrtTwo(svPriTwo.size());
    
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
    
    // Rcpp::print(Rcpp::wrap(myRange));
    
    double wholeTest = (double) myRange / segSize;
    
    if (wholeTest != static_cast<int>(wholeTest))
        ++numCacheSegs;
    
    // Rcpp::print(Rcpp::wrap(numCacheSegs));
    // Rcpp::print(Rcpp::wrap(wholeTest));
    
    int_fast64_t remPrime, timesTwo, maxIndex = myRange + 1;
    bool bKeepGoing;
    
    // Keeps track of which primes will be used in each interval
    // std::deque<std::vector<unsigned long int> > myBuckets(numCacheSegs, 
    //                                                       std::vector<unsigned long int>());

    std::vector<std::vector<unsigned long int>> myBuckets(numCacheSegs,
                                                          std::vector<unsigned long int>());

    for (std::size_t i = 0; i < svPriTwo.size(); ++i) {
        remTest = (lowerBnd % svPriTwo[i]);

        if (remTest == 0) {
            myIndex = svPriTwo[i];
        } else {
            myIndex = (svPriTwo[i] - remTest);
            if ((myIndex % 2) == 0) {myIndex += svPriTwo[i];}
        }

        remTest = (myIndex % sz210) - 1;
        timesTwo = (2 * svPriTwo[i]);
        remPrime = (timesTwo % sz210);
        bKeepGoing = (myIndex < maxIndex);

        while (check210[remTest] && bKeepGoing) {
            myIndex += timesTwo;
            bKeepGoing = (myIndex < maxIndex);
            remTest = remainder210[remTest + remPrime];
        }

        if (bKeepGoing) {
            divTest = (myIndex / segSize);
            myIndex -= (divTest * segSize);
            myBuckets[divTest].push_back(i);
            nextStrtTwo[i] = myIndex;
        }
    }

    int_fast64_t flrMaxNum = segSize * std::floor((double) maxNum / segSize);
    int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
    int_fast64_t myNum = 1 + lowerBnd;

    // vector used for sieving
    std::vector<bool> sieve(segSize, true);
    std::size_t strt = 0;
    Rcpp::print(Rcpp::wrap(svPriOneSize));
    Rcpp::print(Rcpp::wrap(svPriTwo.size()));

    if (minNum != lowerBnd) {
        Rcpp::print(Rcpp::wrap("firstHit"));
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;

            nextStrtOne[i] = (j - segSize);
        }

        for (const auto& elem: myBuckets[0]) {
            myIndex = nextStrtTwo[elem];
            sieve[myIndex] = false;

            // Find the next number divisible by sieve2dPri[i]
            timesTwo = (svPriTwo[elem] * 2);
            myIndex += timesTwo;

            remTest = (myIndex % sz210) - 1;
            remPrime = (timesTwo % sz210);
            bKeepGoing = (myIndex < maxIndex);

            while (check210[remTest] && bKeepGoing) {
                myIndex += timesTwo;
                bKeepGoing = (myIndex < maxIndex);
                remTest = remainder210[remTest + remPrime];
            }

            if (bKeepGoing) {
                divTest = (myIndex / segSize);
                myBuckets[divTest].push_back(elem);
                myIndex -= (divTest * segSize);
                nextStrtTwo[elem] = myIndex;
            }
        }

        if (upperBnd < flrMaxNum) {
            for (std::size_t q = 0; q < numWheelSegs; ++q)
                for (std::size_t w = 0; w < wheelSize; myNum += wheel210[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back((double) myNum);
        } else {
            for (std::size_t q = 0; q < numWheelSegs && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; myNum += wheel210[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back((double) myNum);
        }

        std::fill(sieve.begin(), sieve.end(), true);
        lowerBnd += segSize;
        maxIndex -= segSize;
        ++strt;
    }
    
    Rcpp::print(Rcpp::wrap(flrMaxNum));
    Rcpp::print(Rcpp::wrap(maxNum));

    for (std::size_t v = strt; lowerBnd < flrMaxNum; ++v, lowerBnd += segSize, maxIndex -= segSize) {
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;

            nextStrtOne[i] = (j - segSize);
        }

        for (const auto& elem: myBuckets[v]) {
            myIndex = nextStrtTwo[elem];
            sieve[myIndex] = false;

            // Find the next number divisible by sieve2dPri[i]
            timesTwo = (svPriTwo[elem] * 2);
            myIndex += timesTwo;

            remTest = (myIndex % sz210) - 1;
            remPrime = (timesTwo % sz210);
            bKeepGoing = (myIndex < maxIndex);

            while (check210[remTest] && bKeepGoing) {
                myIndex += timesTwo;
                bKeepGoing = (myIndex < maxIndex);
                remTest = remainder210[remTest + remPrime];
            }

            if (bKeepGoing) {
                divTest = (myIndex / segSize);
                int_fast64_t divPlus = (divTest + v);
                myBuckets[divPlus].push_back(elem);
                myIndex -= (divTest * segSize);
                nextStrtTwo[elem] = myIndex;
            }
        }
        
        for (std::size_t q = 0; q < numWheelSegs; ++q)
            for (std::size_t w = 0; w < wheelSize; myNum += wheel210[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back((double) myNum);

        std::fill(sieve.begin(), sieve.end(), true);
    }
    
    // Get remaining primes that are greater than flrMaxNum and less than maxNum
    if (lowerBnd < maxNum) {
        Rcpp::print(Rcpp::wrap("secondHit"));
        for (std::size_t i = 3; i < svPriOneSize; ++i) {
            int_fast64_t j = nextStrtOne[i];
            for (int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                sieve[j] = false;
        }

        for (const auto& elem: myBuckets[numCacheSegs - 1]) {
            myIndex = nextStrtTwo[elem];
            sieve[myIndex] = false;
        }

        for (std::size_t q = 0; (myNum <= maxNum) && (q < numWheelSegs); ++q)
            for (std::size_t w = 0; (myNum <= maxNum) && (w < wheelSize); myNum += wheel210[w], ++w)
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back((double) myNum);
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

    if (sqrtBound < smallPrimeBase[smlPriBsSize - 1]) {
        if (bAddZero) sievePrimes.push_back(0);
        unsigned long int ind = (bAddTwo) ? 0 : 1;
        
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            sievePrimes.push_back((int_fast32_t) smallPrimeBase[ind]);
        
        if (bAddExtraPrime)
            sievePrimes.push_back((int_fast32_t) smallPrimeBase[ind]);
    } else {
        int sqrtSqrtBound = (int) std::sqrt(sqrtBound);
        std::vector<int_fast64_t> sqrtSievePrimes;
        sqrtSmallPrimes(sqrtSqrtBound, sqrtSievePrimes);
        
        // The number, 225, comes from the observation that the largest prime
        // gap less than 100 million is 219 @ 47,326,693. This is important
        // because we need to guarantee that we obtain the smallest prime
        // greater than sqrt(n). We know that the largest number that this
        // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266
        
        int_fast32_t myLower = 3, myUpper = sqrtBound;
        if (bAddExtraPrime) {myUpper += 225;}
        std::size_t sqrtReserve = EstimatePiPrime(1.0, (double) myUpper);
        sievePrimes.reserve(sqrtReserve);
        
        if (bAddZero) {sievePrimes.push_back(0);}
        if (bAddTwo) {myLower = 1;}
        PrimeSieveMedium(myLower, myUpper, sqrtSievePrimes, sqrtSqrtBound, sievePrimes, 1);
    }
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp2 (SEXP Rb1, SEXP Rb2, SEXP RNumThreads, int nCacheL1, bool bigP, bool medP) {
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
    
    int numThreads, sqrtBound = (int) std::sqrt(myMax);
    int_fast64_t myRange = (int_fast64_t) (myMax - myMin);
    int totalThreads = (int) std::thread::hardware_concurrency();
    
    if ((myRange < 200000) || (totalThreads < 2)) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        Parallel = true;
        CleanConvert::convertPrimitive(RNumThreads, numThreads, "nThreads must be of type numeric or integer");
        if (numThreads > totalThreads) numThreads = totalThreads;
        if (numThreads < 2) Parallel = false;
    }
    
    std::vector<int_fast64_t> sievePrimes;
    double smallLimit = smallPrimeBase[smlPriBsSize - 1];
    smallLimit *= smallLimit;

    // if (myMax >= smallLimit) {
    //     sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
    //     std::vector<int_fast64_t> primes, sievePrimesOne, sievePrimesTwo;
    //     std::size_t ind = 0, limitOne = numSegs * segmentSize;
    //     
    //     // Get the primes that are guaranteed to mark an
    //     // index in the the every segment interval
    //     for ( ; (2 * sievePrimes[ind]) < limitOne; ++ind)
    //         sievePrimesOne.push_back(sievePrimes[ind]);
    //     
    //     // Get the rest
    //     for ( ; ind < sievePrimes.size(); ++ind)
    //         sievePrimesTwo.push_back(sievePrimes[ind]);
    //     
    // 
    //     if (Parallel) {
    //         std::vector<std::thread> myThreads;
    //         std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
    //         int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
    //         int_fast64_t upperBnd = myMin + stepSize - 1;
    // 
    //         for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
    //             myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                    upperBnd, std::ref(sievePrimes),
    //                                    sqrtBound, std::ref(primeList[j]));
    //         }
    // 
    //         myThreads.emplace_back(PrimeSieveBig, lowerBnd,
    //                                (int_fast64_t) myMax, std::ref(sievePrimes),
    //                                sqrtBound, std::ref(primeList[numThreads - 1]));
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
    //         if (numPrimes < INT_MAX) {
    //             Rcpp::NumericVector primes(numPrimes);
    // 
    //             for (int i = 0; i < numThreads; ++i)
    //                 for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
    //                     primes[count] = primeList[i][j];
    // 
    //             return primes;
    //         } else {
    //             return Rcpp::wrap(primeList);
    //         }
    //     } else {
    //         PrimeSieveBig((int_fast64_t) myMin, (int_fast64_t) myMax, sievePrimesOne, sievePrimesTwo, primes);
    //         return Rcpp::wrap(primes);
    //     }
    // }
    //     std::vector<int_fast16_t> small16Primes = sqrt16BasePrimes(sqrtBound);
    //     std::vector<std::vector<int_fast32_t> > primeList(numThreads, std::vector<int_fast32_t>());
    //     int_fast32_t lowerBnd = myMin, stepSize = myRange / numThreads;
    //     int_fast32_t upperBnd = myMin + stepSize - 1;
    // 
        // for (int j = 0; j < (numThreads - 1); ++j) {
        //     myThreads.emplace_back(PrimeSieveSmall, lowerBnd, upperBnd,
        //                            std::ref(small16Primes), sqrtBound, std::ref(primeList[j]));
        //     lowerBnd += stepSize;
        //     upperBnd += stepSize;
        // }
        // 
        // myThreads.emplace_back(PrimeSieveSmall, lowerBnd, (int_fast32_t) myMax,
        //                        std::ref(small16Primes), sqrtBound, std::ref(primeList[numThreads - 1]));
        // 
        // for (auto& thr: myThreads)
        //     thr.join();
        // 
        // for (int i = 0; i < numThreads; ++i) {
        //     numPrimes += primeList[i].size();
        //     sectionSize[i] = primeList[i].size();
        // }
        // 
        // Rcpp::IntegerVector primes(numPrimes);
        // 
        // for (int i = 0; i < numThreads; ++i)
        //     for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
        //         primes[count] = primeList[i][j];
        // 
        // return primes;
    // }
    
    if (bigP) {
        Rcpp::print(Rcpp::wrap("brah"));
        sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
        std::vector<int_fast64_t> sievePrimesOne, sievePrimesTwo;
        std::size_t ind = 0, limitOne = (nCacheL1 * segmentSize);
        
        // Get the primes that are guaranteed to mark an
        // index in the the every segment interval
        for (; (2 * sievePrimes[ind]) < limitOne; ++ind)
            sievePrimesOne.push_back(sievePrimes[ind]);
        
        // Get the rest
        for (; ind < sievePrimes.size(); ++ind)
            sievePrimesTwo.push_back(sievePrimes[ind]);
        
        if (Parallel) {
            if (myRange > 201277440) {
                myRange = 201277440;
            }
            // 
            int_fast64_t tempMax = myMin + myRange;
            std::vector<Rcpp::NumericVector> primeTemp;
            int_fast64_t stepSize = myRange / numThreads;
            int_fast64_t lowerBnd = myMin;
            int_fast64_t upperBnd = lowerBnd + stepSize - 1;
            
            Rcpp::print(Rcpp::wrap(tempMax));
            Rcpp::print(Rcpp::wrap(lowerBnd));
            Rcpp::print(Rcpp::wrap(upperBnd));
            
            Rcpp::print(Rcpp::wrap(myMax));
            
            for (; tempMax <= (int_fast64_t) myMax; tempMax += myRange) {
                // Rcpp::print(Rcpp::wrap(tempMax));
                // Rcpp::print(Rcpp::wrap(lowerBnd));
                std::vector<std::thread> myThreads;
                std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
            
                for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
                    myThreads.emplace_back(PrimeSieveBig, lowerBnd,
                                           upperBnd, std::ref(sievePrimesOne),
                                           std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
                }

                myThreads.emplace_back(PrimeSieveBig, lowerBnd,
                                       (int_fast64_t) tempMax, std::ref(sievePrimesOne),
                                       std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);

                for (auto& thr: myThreads)
                    thr.join();

                std::size_t numPrimes = 0, count = 0;
                std::vector<unsigned long int> sectionSize(numThreads);

                for (int i = 0; i < numThreads; ++i) {
                    numPrimes += primeList[i].size();
                    sectionSize[i] = primeList[i].size();
                }
                
                Rcpp::NumericVector primes(numPrimes);

                for (int i = 0; i < numThreads; ++i)
                    for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                        primes[count] = primeList[i][j];

                primeTemp.push_back(primes);

                lowerBnd = tempMax;
                upperBnd = lowerBnd + stepSize - 1;
            }

            Rcpp::print(Rcpp::wrap("lowerBnd"));
            Rcpp::print(Rcpp::wrap(lowerBnd));
            Rcpp::print(Rcpp::wrap(tempMax));

            stepSize = (myMax - (tempMax - myRange)) / numThreads;
            upperBnd = lowerBnd + stepSize - 1;
            std::vector<std::thread> myThreads;
            std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());

            Rcpp::print(Rcpp::wrap(upperBnd));
            Rcpp::print(Rcpp::wrap(stepSize));

            for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
                myThreads.emplace_back(PrimeSieveBig, lowerBnd,
                                       upperBnd, std::ref(sievePrimesOne),
                                       std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
            }

            myThreads.emplace_back(PrimeSieveBig, lowerBnd,
                                   (int_fast64_t) myMax, std::ref(sievePrimesOne),
                                   std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);

            for (auto& thr: myThreads)
                thr.join();

            std::size_t numPrimes = 0, count = 0;
            std::vector<unsigned long int> sectionSize(numThreads);

            for (int i = 0; i < numThreads; ++i) {
                numPrimes += primeList[i].size();
                sectionSize[i] = primeList[i].size();
            }

            Rcpp::NumericVector primes(numPrimes);

            for (int i = 0; i < numThreads; ++i)
                for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                    primes[count] = primeList[i][j];

            primeTemp.push_back(primes);
            
            return Rcpp::wrap(primeTemp);
            
            // int_fast64_t tempMax = myMin + myRange;
            // std::vector<Rcpp::NumericVector> primeTemp;
            // 
            // 
            // Rcpp::print(Rcpp::wrap("brah2"));
            // std::vector<std::thread> myThreads;
            // std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
            // int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
            // int_fast64_t upperBnd = lowerBnd + stepSize - 1;
            // 
            // for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
            //     myThreads.emplace_back(PrimeSieveBig, lowerBnd,
            //                            upperBnd, std::ref(sievePrimesOne),
            //                            std::ref(sievePrimesTwo), std::ref(primeList[j]), nCacheL1);
            // }
            // 
            // myThreads.emplace_back(PrimeSieveBig, lowerBnd,
            //                        (int_fast64_t) myMax, std::ref(sievePrimesOne),
            //                        std::ref(sievePrimesTwo), std::ref(primeList[numThreads - 1]), nCacheL1);
            // 
            // for (auto& thr: myThreads)
            //     thr.join();
            // 
            // std::size_t numPrimes = 0, count = 0;
            // std::vector<unsigned long int> sectionSize(numThreads);
            // 
            // for (int i = 0; i < numThreads; ++i) {
            //     numPrimes += primeList[i].size();
            //     sectionSize[i] = primeList[i].size();
            // }
            // 
            // if (numPrimes < INT_MAX) {
            //     Rcpp::NumericVector primes(numPrimes);
            // 
            //     for (int i = 0; i < numThreads; ++i)
            //         for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
            //             primes[count] = primeList[i][j];
            // 
            //     return primes;
            // } else {
            //     return Rcpp::wrap(primeList);
            // }
        } else {
            Rcpp::print(Rcpp::wrap("brah3"));
            std::vector<double> primes;
            PrimeSieveBig((int_fast64_t) myMin, (int_fast64_t) myMax,
                          sievePrimesOne, sievePrimesTwo, primes, nCacheL1);
            return Rcpp::wrap(primes);
        }
    } else if (medP) {
        Rcpp::print(Rcpp::wrap("brah20"));
        sqrtBigPrimes(sqrtBound, false, true, false, sievePrimes);
        
        if (Parallel) {
            Rcpp::print(Rcpp::wrap("brah22"));
            
            if (myRange > 1000000000) {
                myRange = 1000000000;
            }
            
            int_fast64_t tempMax = myMin + myRange;
            std::vector<Rcpp::NumericVector> primeTemp;
            int_fast64_t stepSize = myRange / numThreads;
            int_fast64_t lowerBnd = myMin;
            int_fast64_t upperBnd = lowerBnd + stepSize - 1;
            
            for (; tempMax < (int_fast64_t) myMax; tempMax += myRange) {
                Rcpp::print(Rcpp::wrap(tempMax));
                Rcpp::print(Rcpp::wrap(lowerBnd));
                std::vector<std::thread> myThreads;
                std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
                
                for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
                    myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
                                           upperBnd, std::ref(sievePrimes),
                                           sqrtBound, std::ref(primeList[j]), nCacheL1);
                }
                
                myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
                                       (int_fast64_t) tempMax, std::ref(sievePrimes),
                                       sqrtBound, std::ref(primeList[numThreads - 1]), nCacheL1);
                
                for (auto& thr: myThreads)
                    thr.join();
                
                std::size_t numPrimes = 0, count = 0;
                std::vector<unsigned long int> sectionSize(numThreads);
                
                for (int i = 0; i < numThreads; ++i) {
                    numPrimes += primeList[i].size();
                    sectionSize[i] = primeList[i].size();
                }
                
                Rcpp::NumericVector primes(numPrimes);
                
                for (int i = 0; i < numThreads; ++i)
                    for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                        primes[count] = primeList[i][j];
                
                primeTemp.push_back(primes);
                
                lowerBnd = tempMax;
                upperBnd = lowerBnd + stepSize - 1;
            }
            
            Rcpp::print(Rcpp::wrap("lowerBnd"));
            Rcpp::print(Rcpp::wrap(lowerBnd));
            Rcpp::print(Rcpp::wrap(tempMax));
            
            stepSize = (myMax - (tempMax - myRange)) / numThreads;
            upperBnd = lowerBnd + stepSize - 1;
            std::vector<std::thread> myThreads;
            std::vector<std::vector<double>> primeList(numThreads, std::vector<double>());
            
            Rcpp::print(Rcpp::wrap(upperBnd));
            Rcpp::print(Rcpp::wrap(stepSize));

            for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
                myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
                                       upperBnd, std::ref(sievePrimes),
                                       sqrtBound, std::ref(primeList[j]), nCacheL1);
            }

            myThreads.emplace_back(PrimeSieveMedium<double>, lowerBnd,
                                   (int_fast64_t) myMax, std::ref(sievePrimes),
                                   sqrtBound, std::ref(primeList[numThreads - 1]), nCacheL1);

            for (auto& thr: myThreads)
                thr.join();

            std::size_t numPrimes = 0, count = 0;
            std::vector<unsigned long int> sectionSize(numThreads);

            for (int i = 0; i < numThreads; ++i) {
                numPrimes += primeList[i].size();
                sectionSize[i] = primeList[i].size();
            }

            Rcpp::NumericVector primes(numPrimes);

            for (int i = 0; i < numThreads; ++i)
                for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                    primes[count] = primeList[i][j];

            primeTemp.push_back(primes);
            
            return Rcpp::wrap(primeTemp);
            // if (numPrimes < INT_MAX) {
                // Rcpp::NumericVector primes(numPrimes);
                // 
                // for (int i = 0; i < numThreads; ++i)
                //     for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                //         primes[count] = primeList[i][j];
                // 
                // return primes;
            // } else {
            //     return Rcpp::wrap(primeList);
            // }
        } else {
            Rcpp::print(Rcpp::wrap("brah23"));
            std::vector<double> primes;
            PrimeSieveMedium((int_fast64_t) myMin, (int_fast64_t) myMax,
                          sievePrimes, sqrtBound, primes, nCacheL1);
            return Rcpp::wrap(primes);
        }
    }
    
    sqrtSmallPrimes(sqrtBound, sievePrimes);
    std::vector<int_fast32_t> primes;
    
    if (Parallel) {
        std::vector<std::thread> myThreads;
        std::vector<std::vector<int_fast32_t> > primeList(numThreads, std::vector<int_fast32_t>());
        int_fast64_t lowerBnd = myMin, stepSize = myRange / numThreads;
        int_fast64_t upperBnd = myMin + stepSize - 1;
        
        for (int j = 0; j < (numThreads - 1); ++j, lowerBnd += stepSize, upperBnd += stepSize) {
            myThreads.emplace_back(PrimeSieveSmall, lowerBnd, upperBnd,
                                   std::ref(sievePrimes), sqrtBound, std::ref(primeList[j]));
        }
        
        myThreads.emplace_back(PrimeSieveSmall, lowerBnd, (int_fast64_t) myMax,
                               std::ref(sievePrimes), sqrtBound, std::ref(primeList[numThreads - 1]));
        
        for (auto& thr: myThreads)
            thr.join();
        
        std::size_t numPrimes = 0, count = 0;
        std::vector<unsigned long int> sectionSize(numThreads);
        
        for (int i = 0; i < numThreads; ++i) {
            numPrimes += primeList[i].size();
            sectionSize[i] = primeList[i].size();
        }
        
        Rcpp::IntegerVector primes(numPrimes);
        
        for (int i = 0; i < numThreads; ++i)
            for (std::size_t j = 0; j < sectionSize[i]; ++j, ++count)
                primes[count] = primeList[i][j];
        
        return primes;
    } else {
        PrimeSieveSmall((int_fast64_t) myMin, (int_fast64_t) myMax, sievePrimes, sqrtBound, primes);
        return Rcpp::wrap(primes);
    }
}
