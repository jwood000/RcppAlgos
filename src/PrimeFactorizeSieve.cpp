#include "NumbersUtils/MotleyStartIndex.h"
#include "NumbersUtils/libdivide.h"
#include <vector>
#include <cmath>

namespace MotleyPrimes {

    template <typename T>
    void PrimeFactorizationSieve(T m, T retN, T offsetStrt,
                                 const std::vector<T> &primes,
                                 std::vector<std::vector<T>> &primeList) {

        const T n = retN;
        const T myRange = (n - m) + 1;
        T myNum = m;

        const double myLogN = std::log(static_cast<double>(n));

        if (n > 3) {
            std::vector<std::uint8_t> myMemory(myRange, 1u);
            const T sqrtBound = std::sqrt(static_cast<double>(retN));

            for (auto p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                const std::size_t limit = myLogN / std::log(static_cast<double>(*p));

                if (m < 2) {
                    for (std::size_t i = 1; i <= limit; ++i) {
                        const T myStep = static_cast<T>(std::pow(*p, i));

                        for (T j = (myStep - 1); j < n; j += myStep) {
                            ++myMemory[j];
                        }
                    }
                } else {
                    for (std::size_t i = 1; i <= limit; ++i) {
                        const T myStep = static_cast<T>(std::pow(*p, i));

                        for (T j = getStartIndexPowP(m, myStep, *p);
                             j < myRange; j += myStep) {
                            ++myMemory[j];
                        }
                    }
                }
            }

            std::vector<uint8_t>::iterator myMalloc;
            typename std::vector<std::vector<T>>::iterator it2d;
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
                for (auto p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                    const std::size_t limit = myLogN / std::log(static_cast<double>(*p));
                    const libdivide::divider<T> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        const T myStep = static_cast<T>(std::pow(*p, i));

                        for (T j = (myStep - 1); j < n; j += myStep) {
                            if (primeList[j].back() > *p) {
                                primeList[j].back() /= fastDiv;
                                primeList[j].insert(primeList[j].end() - 1, *p);
                            }
                        }
                    }
                }
            } else {
                const T offsetRange = myRange + offsetStrt;

                for (auto p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                    const std::size_t limit = myLogN / std::log(static_cast<double>(*p));
                    const libdivide::divider<T> fastDiv(*p);

                    for (std::size_t i = 1; i <= limit; ++i) {
                        const T myStep = static_cast<T>(std::pow(*p, i));
                        const T myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);

                        for (T j = myStart; j < offsetRange; j += myStep) {
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

            for (int i = strt; i < myRange; ++i, ++myNum) {
                primeList[i].push_back(myNum);
            }
        }
    }
}

template void MotleyPrimes::PrimeFactorizationSieve(
        std::int64_t, std::int64_t, std::int64_t,
        const std::vector<std::int64_t>&,
        std::vector<std::vector<std::int64_t>>&
);

template void MotleyPrimes::PrimeFactorizationSieve(
        int, int, int, const std::vector<int>&,
        std::vector<std::vector<int>>&
);
