#include "MotleyStartIndex.h"
#include <libdivide.h>
#include <cmath>

namespace MotleyPrimes {

    template <typename typeInt>
    void PrimeFactorizationSieve(typeInt m, typeInt retN, typeInt offsetStrt,
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
                const std::size_t limit = static_cast<std::size_t>(trunc(myLogN / std::log(*p)));
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
            const auto itEnd = primeList.begin() + offsetStrt + myRange;
            
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
                    const std::size_t limit = static_cast<std::size_t>(trunc(myLogN / std::log(*p)));
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
                    const std::size_t limit = static_cast<std::size_t>(trunc(myLogN / std::log(*p)));
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
}


template void MotleyPrimes::PrimeFactorizationSieve(std::int64_t, std::int64_t, std::int64_t,
                                                    const std::vector<std::int64_t>&,
                                                    std::vector<std::vector<std::int64_t>>&);

template void MotleyPrimes::PrimeFactorizationSieve(int, int, int,
                                                    const std::vector<int>&,
                                                    std::vector<std::vector<int>>&);
