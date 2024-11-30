#include "CppConvert.h"

#include "NumbersUtils/PrimeFactorizeSieve.h"
#include "NumbersUtils/EulerPhiSieve.h"
#include "NumbersUtils/Eratosthenes.h"
#include <type_traits>
#include <thread>

template <typename T, typename U>
void MotleyMain(T myMin, U myMax, bool IsEuler,
                U* EulerPhis, std::vector<T> &numSeq,
                std::vector<std::vector<T>> &primeList,
                int nThreads, int maxThreads) {

    bool Parallel = false;
    std::int_fast64_t myRange = (myMax - myMin) + 1;
    T offsetStrt = 0;

    if (nThreads > 1 && maxThreads > 1 && myRange >= 20000) {
        Parallel = true;
        if (nThreads > maxThreads) nThreads = maxThreads;

        // Ensure that each thread has at least 10000
        if ((myRange / nThreads) < 10000) nThreads = myRange / 10000;
    }

    std::vector<T> primes;
    int sqrtBound = std::sqrt(static_cast<double>(myMax));
    PrimeSieve::sqrtBigPrimes(sqrtBound, false, true, true, primes);

    if (Parallel) {
        std::vector<std::thread> threads;
        T lower = myMin;
        T chunkSize = myRange / nThreads;
        U upper = lower + chunkSize - 1;

        for (int ind = 0; ind < (nThreads - 1); offsetStrt += chunkSize,
             lower = (upper + 1), upper += chunkSize, ++ind) {
            if (IsEuler) {
                threads.emplace_back(std::cref(MotleyPrimes::EulerPhiSieve<T, U>),
                                     lower, upper, offsetStrt, std::ref(primes),
                                     std::ref(numSeq), EulerPhis);
            } else {
                threads.emplace_back(std::cref(MotleyPrimes::PrimeFactorizationSieve<T>),
                                     lower, static_cast<T>(upper), offsetStrt,
                                     std::cref(primes), std::ref(primeList));
            }
        }

        if (IsEuler) {
            threads.emplace_back(std::cref(MotleyPrimes::EulerPhiSieve<T, U>),
                                 lower, myMax, offsetStrt, std::ref(primes),
                                 std::ref(numSeq), EulerPhis);
        } else {
            threads.emplace_back(std::cref(MotleyPrimes::PrimeFactorizationSieve<T>),
                                 lower, static_cast<T>(myMax), offsetStrt,
                                 std::cref(primes), std::ref(primeList));
        }

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        if (IsEuler) {
            MotleyPrimes::EulerPhiSieve(myMin, myMax, offsetStrt,
                                        primes, numSeq, EulerPhis);
        } else {
            MotleyPrimes::PrimeFactorizationSieve(myMin, static_cast<T>(myMax),
                                                  offsetStrt, primes, primeList);
        }
    }
}

template <typename T, typename U>
SEXP GlueMotley(T myMin, U myMax, bool IsEuler,
                bool keepNames, int nThreads, int maxThreads) {

    std::size_t myRange = (myMax - myMin) + 1;

    if (IsEuler) {
        if (std::is_integral<U>::value) {
            std::vector<std::vector<int>> tempList;
            std::vector<int> numSeq(myRange);
            cpp11::r_vector<int> EulerPhis(Rf_allocVector(INTSXP, myRange));
            int* ptrEuler = INTEGER(EulerPhis);

            MotleyMain(
                static_cast<int>(myMin), static_cast<int>(myMax),
                IsEuler, ptrEuler, numSeq, tempList, nThreads, maxThreads
            );

            if (keepNames) {
                CppConvert::SetNames(
                    EulerPhis, static_cast<int>(myMin), static_cast<int>(myMax)
                );
            }

            return EulerPhis;
        } else {
            std::vector<std::vector<std::int64_t>> tempList;
            std::vector<std::int64_t> numSeq(myRange);
            cpp11::r_vector<double> EulerPhis(Rf_allocVector(REALSXP, myRange));
            double* ptrEuler = REAL(EulerPhis);

            MotleyMain(
                static_cast<std::int64_t>(myMin), static_cast<double>(myMax),
                IsEuler, ptrEuler, numSeq, tempList, nThreads, maxThreads
            );

            if (keepNames) {
                CppConvert::SetNames(
                    EulerPhis, static_cast<double>(myMin),
                    static_cast<double>(myMax)
                );
            }

            return EulerPhis;
        }
    } else {
        std::vector<std::vector<T>>
            primeList(myRange, std::vector<T>());
        U* tempEuler = nullptr;
        std::vector<T> tempVec;

        MotleyMain(myMin, myMax, IsEuler, tempEuler,
                   tempVec, primeList, nThreads, maxThreads);
        cpp11::writable::list myList(myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            myList[i] = cpp11::writable::r_vector<U>(primeList[i]);
        }

        if (keepNames) {
            CppConvert::SetNames(
                myList, static_cast<U>(myMin), static_cast<U>(myMax)
            );
        }
        return myList;
    }
}

[[cpp11::register]]
SEXP MotleyContainer(SEXP Rb1, SEXP Rb2, SEXP RIsEuler, SEXP RNamed,
                     SEXP RNumThreads, SEXP RmaxThreads) {

    double bound1 = 0;
    double bound2 = 0;

    double myMin;
    double myMax;

    int nThreads = 1;
    int maxThreads = 1;

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");
    const bool IsEuler = CppConvert::convertFlag(RIsEuler, "IsEuler");
    const std::string namedObject = (IsEuler) ? "namedVector" : "namedList";
    bool IsNamed = CppConvert::convertFlag(RNamed, namedObject);
    CppConvert::convertPrimitive(Rb1, bound1, VecType::Numeric, "bound1");

    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CppConvert::convertPrimitive(Rb2, bound2, VecType::Numeric, "bound2");
    }

    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }

    if (myMax < 2) {
        using namespace cpp11::literals;
        if (IsEuler) {
            cpp11::writable::r_vector<int> res({1});
            if (IsNamed) Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            return res;
        } else {
            cpp11::writable::list res({
                cpp11::r_vector<int>(Rf_allocVector(INTSXP, 0))
            });
            if (IsNamed) Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            return res;
        }
    }

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    if (myMax > std::numeric_limits<int>::max()) {
        return GlueMotley(static_cast<std::int64_t>(myMin), myMax,
                          IsEuler, IsNamed, nThreads, maxThreads);
    } else {
        return GlueMotley(static_cast<int>(myMin), static_cast<int>(myMax),
                          IsEuler, IsNamed, nThreads, maxThreads);
    }
}
