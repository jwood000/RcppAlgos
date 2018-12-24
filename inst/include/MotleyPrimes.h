#ifndef MOTLEY_PRIMES_H
#define MOTLEY_PRIMES_H

#include "EratosthenesSieve.h"
#include <libdivide.h>

// This file contains functions for generating the prime factorization
// and applying the euler phi function over a range of numbers.

namespace MotleyPrimes {

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
    void PrimeFactorizationSieve(typeInt m, typeReturn retN, typeInt offsetStrt,
                                 const std::vector<typeInt> &primes,
                                 std::vector<std::vector<typeReturn>> &primeList) {
        
        typeInt n = (typeInt) retN;
        typeInt myRange = n;
        myRange += (1 - m);
        
        typeInt myStep, myStart, myNum = m;
        double myLogN = std::log(n);
        unsigned long int limit;
        
        if (n > 3) {
            typename std::vector<typeInt>::const_iterator p;
            std::vector<uint8_t> myMemory(myRange, 1u);
            typeInt sqrtBound = static_cast<typeInt>(sqrt(retN));
            
            for (p = primes.begin(); (*p) <= sqrtBound; ++p) {
                limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
                if (m < 2) {
                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));
                        for (typeInt j = (myStep - 1); j < n; j += myStep)
                            ++myMemory[j];
                    }
                } else {
                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));
                        myStart = getStartIndexPowP(m, myStep, *p);
                        for (typeInt j = myStart; j < myRange; j += myStep)
                            ++myMemory[j];
                    }
                }
            }
            
            std::vector<uint8_t>::iterator myMalloc;
            typename std::vector<std::vector<typeReturn>>::iterator it2d, itEnd;
            itEnd = primeList.begin() + offsetStrt + myRange;
            
            if (myNum < 2) {
                ++myNum;
                it2d = primeList.begin() + 1;
                myMalloc = myMemory.begin() + 1;
            } else {
                it2d = primeList.begin() + offsetStrt;
                myMalloc = myMemory.begin();
            }

            for (; it2d != itEnd; ++it2d, ++myNum, ++myMalloc) {
                it2d->reserve(*myMalloc);
                it2d->push_back(static_cast<typeReturn>(myNum));
            }
            
            if (m < 2) {
                for (p = primes.begin(); (*p) <= sqrtBound; ++p) {
                    limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
                    libdivide::divider<typeInt> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));

                        for (typeInt j = (myStep - 1); j < n; j += myStep) {
                            if (primeList[j].back() > *p) {
                                myNum = static_cast<typeInt>(primeList[j].back());
                                myNum /= fastDiv;
                                primeList[j].back() = static_cast<typeReturn>(myNum);
                                primeList[j].insert(primeList[j].end() - 1, static_cast<typeReturn>(*p));
                            }
                        }
                    }
                }
            } else {
                typeInt offsetRange = myRange + offsetStrt;
                
                for (p = primes.begin(); (*p) <= sqrtBound; ++p) {
                    limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
                    libdivide::divider<typeInt> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));
                        myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);

                        for (typeInt j = myStart; j < offsetRange; j += myStep) {
                            if (primeList[j].back() > *p) {
                                myNum = static_cast<typeInt>(primeList[j].back());
                                myNum /= fastDiv;
                                primeList[j].back() = static_cast<typeReturn>(myNum);
                                primeList[j].insert(primeList[j].end() - 1, static_cast<typeReturn>(*p));
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
                primeList[i].push_back(static_cast<typeReturn>(myNum));
        }
    }
    
    template <typename typeInt, typename typeReturn, typename typeRcpp>
    void EulerPhiSieve(typeInt m, typeReturn retN, typeInt offsetStrt,
                       const std::vector<typeInt> &primes,
                       std::vector<typeInt> &numSeq,
                       typeRcpp &EulerPhis) {
        
        typeInt n = static_cast<typeInt>(retN);
        typeInt myRange = n;
        myRange += (1 - m);
        
        typeInt myNum = m;
        double myLogN = std::log(n);
        unsigned long int limit;
        typeReturn retNum = static_cast<typeReturn>(m);
        
        for (std::size_t i = offsetStrt; retNum <= retN; ++retNum, ++i) {
            EulerPhis[i] = retNum;
            numSeq[i] = static_cast<typeInt>(retNum);
        }
        
        if (m < 2) {
            bool tempPar = false;
            std::vector<typeInt> fullPrimes;
            int_fast64_t intMin = static_cast<int_fast64_t>(m);
            int_fast64_t intMax = static_cast<int_fast64_t>(retN);
            std::vector<std::vector<typeInt>> tempList;
            PrimeSieve::PrimeSieveMaster(intMin, intMax, fullPrimes, tempList, tempPar);
            typename std::vector<typeInt>::iterator p;
            
            for (p = fullPrimes.begin(); p < fullPrimes.end(); ++p) {
                libdivide::divider<typeInt> fastDiv(*p);
                for (typeInt j = (*p - 1); j < n; j += *p) {
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
            }
        } else if (n > 3) {
            typename std::vector<typeInt>::const_iterator p;
            typeInt sqrtBound = static_cast<typeInt>(sqrt(retN));
            typeInt myStart, myStep, offsetRange = myRange + offsetStrt;
            
            for (p = primes.begin(); (*p) <= sqrtBound; ++p) {
                limit = static_cast<unsigned long int>(myLogN / std::log(*p));
                myStart = offsetStrt + getStartIndexPowP(m, *p, *p);
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (typeInt j = myStart; j < offsetRange; j += *p) {
                    numSeq[j] /= fastDiv;
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
                
                for (std::size_t i = 2; i <= limit; ++i) {
                    myStep = static_cast<typeInt>(std::pow(*p, i));
                    myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);
                    
                    for (typeInt j = myStart; j < offsetRange; j += myStep)
                        numSeq[j] /= fastDiv;
                }
            }
            
            for (typeInt i = offsetStrt; i < offsetRange; ++i) {
                if (numSeq[i] > 1) {
                    myNum = static_cast<typeInt>(EulerPhis[i]);
                    myNum /= numSeq[i];
                    EulerPhis[i] -= static_cast<typeReturn>(myNum);
                }
            }
            
        } else { // edge case where m,n = 2 or 3
            for (int i = 0; i < myRange; ++i)
                --EulerPhis[i];
        }
    }
    
    template <typename typeInt, typename typeReturn, typename typeRcpp>
    void MotleyMaster(typeInt myMin, typeReturn myMax, bool isEuler,
                      typeRcpp &EulerPhis, std::vector<typeInt> &numSeq,
                      std::vector<std::vector<typeReturn>> &primeList,
                      int nThreads, int maxThreads) {
        
        bool Parallel = false;
        int_fast64_t myRange = (myMax - myMin) + 1;
        typeInt offsetStrt = 0;
        
        if (nThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) {nThreads = maxThreads;}
        }
        
        if ((maxThreads < 2) || (myRange < 10000)) {Parallel = false;}
        int sqrtBound = std::sqrt(myMax);
        std::vector<typeInt> primes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, false, true, true, primes);

        if (Parallel) {
            std::size_t ind = 0u;
            std::vector<std::thread> myThreads;
            typeInt lowerBnd = myMin;
            typeInt chunkSize = myRange / nThreads;
            typeReturn upperBnd = lowerBnd + chunkSize - 1;
            
            for (; ind < (nThreads - 1); offsetStrt += chunkSize, 
                 lowerBnd = (upperBnd + 1), upperBnd += chunkSize, ++ind) {
                if (isEuler) {
                    myThreads.emplace_back(EulerPhiSieve<typeInt, typeReturn, typeRcpp>,
                                           lowerBnd, upperBnd, offsetStrt, 
                                           std::ref(primes), std::ref(numSeq),
                                           std::ref(EulerPhis));
                } else {
                    myThreads.emplace_back(PrimeFactorizationSieve<typeInt, typeReturn>,
                                           lowerBnd, upperBnd, offsetStrt, 
                                           std::ref(primes), std::ref(primeList));
                }
            }
            
            if (isEuler) {
                myThreads.emplace_back(EulerPhiSieve<typeInt, typeReturn, typeRcpp>,
                                       lowerBnd, myMax, offsetStrt,
                                       std::ref(primes), std::ref(numSeq),
                                       std::ref(EulerPhis));
            } else {
                myThreads.emplace_back(PrimeFactorizationSieve<typeInt, typeReturn>,
                                       lowerBnd, myMax, offsetStrt,
                                       std::ref(primes), std::ref(primeList));
            }

            for (auto &thr: myThreads)
                thr.join();
            
        } else {
            if (isEuler) {
                EulerPhiSieve(myMin, myMax, offsetStrt, primes, numSeq, EulerPhis);
            } else {
                PrimeFactorizationSieve(myMin, myMax, offsetStrt, primes, primeList);
            }
        }
    }
}

#endif
