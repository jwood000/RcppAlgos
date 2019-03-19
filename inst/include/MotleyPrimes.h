#ifndef MOTLEY_PRIMES_H
#define MOTLEY_PRIMES_H

#include "Eratosthenes.h"
#include <libdivide.h>

// This file contains functions for generating the prime factorization
// and applying the euler phi function over a range of numbers.

namespace MotleyPrimes {

    // This function is slightly different than the getStartingIndex
    // in the DivisorsContainer.cpp file. The step passed in this
    // function is a power of prime and requires an additional check
    // (i.e. else if (myPrime < lowerB)).
    template <typename typeInt>
    inline typeInt getStartIndexPowP(const typeInt lowerB, const typeInt step, 
                                     const typeInt myPrime) {
        
        typeInt retStrt;
        const typeInt remTest = lowerB % step;
        
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
    void PrimeFactorizationSieve(const typeInt m, const typeInt retN, const typeInt offsetStrt,
                                 const std::vector<typeInt> &primes,
                                 std::vector<std::vector<typeInt>> &primeList) {
        
        const typeInt n = (typeInt) retN;
        const typeInt myRange = (n - m) + 1;
        
        typeInt myStep, myStart, myNum = m;
        const double myLogN = std::log(n);
        
        if (n > 3) {
            typename std::vector<typeInt>::const_iterator p;
            std::vector<uint8_t> myMemory(myRange, 1u);
            const typeInt sqrtBound = static_cast<typeInt>(std::sqrt(static_cast<double>(retN)));
            
            for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                const unsigned long int limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
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
            typename std::vector<std::vector<typeInt>>::iterator it2d;
            const typename std::vector<std::vector<typeInt>>::iterator itEnd = primeList.begin() + offsetStrt + myRange;
            
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
                it2d->push_back(myNum);
            }
            
            if (m < 2) {
                for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                    const unsigned long int limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
                    const libdivide::divider<typeInt> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));

                        for (typeInt j = (myStep - 1); j < n; j += myStep) {
                            if (primeList[j].back() > *p) {
                                primeList[j].back() /= fastDiv;
                                primeList[j].insert(primeList[j].end() - 1, *p);
                            }
                        }
                    }
                }
            } else {
                const typeInt offsetRange = myRange + offsetStrt;
                
                for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                    const unsigned long int limit = static_cast<unsigned long int>(trunc(myLogN / std::log(*p)));
                    const libdivide::divider<typeInt> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        myStep = static_cast<typeInt>(std::pow(*p, i));
                        myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);

                        for (typeInt j = myStart; j < offsetRange; j += myStep) {
                            if (primeList[j].back() > *p) {
                                primeList[j].back() /= fastDiv;
                                primeList[j].insert(primeList[j].end() - 1, *p);
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
                primeList[i].push_back(myNum);
        }
    }
    
    template <typename typeInt, typename typeReturn, typename typeRcpp>
    void EulerPhiSieve(const typeInt m, const typeReturn retN, const typeInt offsetStrt,
                       const std::vector<typeInt> &primes,
                       std::vector<typeInt> &numSeq,
                       typeRcpp &EulerPhis) {
        
        const typeInt n = static_cast<typeInt>(retN);
        const typeInt myRange = (n - m) + 1;
        
        typeInt myNum = m;
        const double myLogN = std::log(n);
        typeReturn retNum = static_cast<typeReturn>(m);
        
        for (std::size_t i = offsetStrt; retNum <= retN; ++retNum, ++i) {
            EulerPhis[i] = retNum;
            numSeq[i] = static_cast<typeInt>(retNum);
        }
        
        if (m < 2) {
            bool tempPar = false;
            std::vector<typeInt> fullPrimes;
            const int_fast64_t intMin = static_cast<int_fast64_t>(m);
            const int_fast64_t intMax = static_cast<int_fast64_t>(retN);
            std::vector<std::vector<typeInt>> tempList;
            PrimeSieve::PrimeSieveMaster(intMin, intMax, fullPrimes, tempList, tempPar);
            typename std::vector<typeInt>::iterator p;
            
            for (p = fullPrimes.begin(); p < fullPrimes.end(); ++p) {
                const libdivide::divider<typeInt> fastDiv(*p);
                for (typeInt j = (*p - 1); j < n; j += *p) {
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
            }
        } else if (n > 3) {
            typename std::vector<typeInt>::const_iterator p;
            const typeInt sqrtBound = static_cast<typeInt>(std::sqrt(static_cast<double>(retN)));
            const typeInt offsetRange = myRange + offsetStrt;
            
            for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                const unsigned long int limit = static_cast<unsigned long int>(myLogN / std::log(*p));
                const typeInt myStart = offsetStrt + getStartIndexPowP(m, *p, *p);
                const libdivide::divider<typeInt> fastDiv(*p);
                
                for (typeInt j = myStart; j < offsetRange; j += *p) {
                    numSeq[j] /= fastDiv;
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
                
                for (std::size_t i = 2; i <= limit; ++i) {
                    const typeInt myStep = static_cast<typeInt>(std::pow(*p, i));
                    const typeInt myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);
                    
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
                      std::vector<std::vector<typeInt>> &primeList,
                      int nThreads, int maxThreads) {
        
        bool Parallel = false;
        int_fast64_t myRange = (myMax - myMin) + 1;
        typeInt offsetStrt = 0;
        
        if (nThreads > 1 && maxThreads > 1 && myRange >= 20000) {
            Parallel = true;
            
            if (nThreads > maxThreads)
                nThreads = maxThreads;
            
            // Ensure that each thread has at least 10000
            if ((myRange / nThreads) < 10000)
                nThreads = myRange / 10000;
        }
        
        int sqrtBound = std::sqrt(static_cast<double>(myMax));
        std::vector<typeInt> primes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, false, true, true, primes);

        if (Parallel) {
            RcppThread::ThreadPool pool(nThreads);
            typeInt lower = myMin;
            typeInt chunkSize = myRange / nThreads;
            typeReturn upper = lower + chunkSize - 1;
            
            for (int ind = 0; ind < (nThreads - 1); offsetStrt += chunkSize, 
                                lower = (upper + 1), upper += chunkSize, ++ind) {
                if (isEuler) {
                    pool.push(std::cref(EulerPhiSieve<typeInt, typeReturn, typeRcpp>), lower,
                              upper, offsetStrt, std::ref(primes), std::ref(numSeq), std::ref(EulerPhis));
                } else {
                    pool.push(std::cref(PrimeFactorizationSieve<typeInt>), 
                              lower, static_cast<typeInt>(upper), offsetStrt,
                              std::ref(primes), std::ref(primeList));
                }
            }
            
            if (isEuler) {
                pool.push(std::cref(EulerPhiSieve<typeInt, typeReturn, typeRcpp>), lower,
                          myMax, offsetStrt, std::ref(primes), std::ref(numSeq), std::ref(EulerPhis));
            } else {
                pool.push(std::cref(PrimeFactorizationSieve<typeInt>), 
                          lower, static_cast<typeInt>(myMax), offsetStrt,
                          std::ref(primes), std::ref(primeList));
            }

            pool.join();
            
        } else {
            if (isEuler)
                EulerPhiSieve(myMin, myMax, offsetStrt, primes, numSeq, EulerPhis);
            else
                PrimeFactorizationSieve(myMin, static_cast<typeInt>(myMax), offsetStrt, primes, primeList);
        }
    }
}

#endif
