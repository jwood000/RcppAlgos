#include "NumbersUtils/MotleyStartIndex.h"
#include "NumbersUtils/Eratosthenes.h"
#include "NumbersUtils/libdivide.h"
#include <cmath>

namespace MotleyPrimes {

    template <typename T, typename U>
    void EulerPhiSieve(T m, U retN, T offsetStrt,
                       const std::vector<T> &primes,
                       std::vector<T> &numSeq,
                       U* EulerPhis) {

        const T n = static_cast<T>(retN);
        const T myRange = (n - m) + 1;

        const double myLogN = std::log(static_cast<double>(n));
        U retNum = static_cast<U>(m);

        for (std::size_t i = offsetStrt; retNum <= retN; ++retNum, ++i) {
            EulerPhis[i] = retNum;
            numSeq[i] = static_cast<T>(retNum);
        }

        if (m < 2) {
            bool tempPar = false;
            std::vector<T> fullPrimes;
            const std::int_fast64_t intMin = static_cast<std::int_fast64_t>(m);
            const std::int_fast64_t intMax = static_cast<std::int_fast64_t>(retN);
            std::vector<std::vector<T>> tempList;
            PrimeSieve::PrimeSieveMain(tempList, fullPrimes, intMin, intMax, tempPar);

            for (auto p: fullPrimes) {
                const libdivide::divider<T> fastDiv(p);
                for (T j = (p - 1); j < n; j += p) {
                    T myNum = static_cast<T>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<U>(myNum);
                }
            }
        } else if (n > 3) {
            typename std::vector<T>::const_iterator p;
            const T sqrtBound = static_cast<T>(std::sqrt(static_cast<double>(retN)));
            const T offsetRange = myRange + offsetStrt;

            for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                const std::size_t limit = static_cast<std::size_t>(myLogN / std::log(*p));
                const T myStart = offsetStrt + getStartIndexPowP(m, *p, *p);
                const libdivide::divider<T> fastDiv(*p);

                for (T j = myStart; j < offsetRange; j += *p) {
                    numSeq[j] /= fastDiv;
                    T myNum = static_cast<T>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<U>(myNum);
                }

                for (std::size_t i = 2; i <= limit; ++i) {
                    const T myStep = static_cast<T>(std::pow(*p, i));
                    const T myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);

                    for (T j = myStart; j < offsetRange; j += myStep) {
                        numSeq[j] /= fastDiv;
                    }
                }
            }

            for (T i = offsetStrt; i < offsetRange; ++i) {
                if (numSeq[i] > 1) {
                    T myNum = static_cast<T>(EulerPhis[i]);
                    myNum /= numSeq[i];
                    EulerPhis[i] -= static_cast<U>(myNum);
                }
            }

        } else { // edge case where m,n = 2 or 3
            for (int i = 0; i < myRange; ++i) {
                --EulerPhis[i];
            }
        }
    }
}

template void MotleyPrimes::EulerPhiSieve(int, int, int,
                                          const std::vector<int>&,
                                          std::vector<int>&, int*);

template void MotleyPrimes::EulerPhiSieve(std::int64_t, double, std::int64_t,
                                          const std::vector<std::int64_t>&,
                                          std::vector<std::int64_t>&, double*);
