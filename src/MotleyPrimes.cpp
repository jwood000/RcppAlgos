#include "NumbersUtils/PrimeFactorizeSieve.h"
#include "NumbersUtils/EulerPhiSieve.h"
#include "NumbersUtils/Eratosthenes.h"
#include "CleanConvert.h"
#include "SetUpUtils.h"
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
            MotleyPrimes::EulerPhiSieve(myMin, myMax, offsetStrt, primes, numSeq, EulerPhis);
        } else {
            MotleyPrimes::PrimeFactorizationSieve(myMin, static_cast<T>(myMax),
                                                  offsetStrt, primes, primeList);
        }
    }
}

SEXP GlueIntMotley(int myMin, int myMax, bool IsEuler,
                   bool keepNames, int nThreads, int maxThreads) {

    int numUnprotects = 1;
    std::size_t myRange = (myMax - myMin) + 1;

    if (IsEuler) {
        std::vector<std::vector<int>> tempList;
        cpp11::sexp EulerPhis = Rf_allocVector(INTSXP, myRange);
        int* ptrEuler  = INTEGER(EulerPhis);
        std::vector<int> numSeq(myRange);

        MotleyMain(myMin, myMax, IsEuler, ptrEuler,
                   numSeq, tempList, nThreads, maxThreads);

        if (keepNames) {
            ++numUnprotects;
            SetIntNames(EulerPhis, myRange, myMin, myMax);
        }

        return EulerPhis;
    } else {
        std::vector<std::vector<int>>
            primeList(myRange, std::vector<int>());

        int* tempEuler = nullptr;
        std::vector<int> tempVec;
        MotleyMain(myMin, myMax, IsEuler, tempEuler,
                   tempVec, primeList, nThreads, maxThreads);

        cpp11::sexp myList = Rf_allocVector(VECSXP, myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            SET_VECTOR_ELT(myList, i, GetIntVec(primeList[i]));
        }

        if (keepNames) {
            ++numUnprotects;
            SetIntNames(myList, myRange, myMin, myMax);
        }

        return myList;
    }
}

SEXP GlueDblMotley(std::int64_t myMin, double myMax, bool IsEuler,
                   bool keepNames, int nThreads, int maxThreads) {

    int numUnprotects = 1;
    std::size_t myRange = (myMax - myMin) + 1;

    if (IsEuler) {
        std::vector<std::vector<std::int64_t>> tempList;
        cpp11::sexp EulerPhis = Rf_allocVector(REALSXP, myRange);
        double* ptrEuler = REAL(EulerPhis);
        std::vector<std::int64_t> numSeq(myRange);

        MotleyMain(myMin, myMax, IsEuler, ptrEuler,
                   numSeq, tempList, nThreads, maxThreads);

        if (keepNames) {
            ++numUnprotects;
            SetDblNames(EulerPhis, myRange,
                        static_cast<double>(myMin), myMax);
        }

        return EulerPhis;
    } else {
        std::vector<std::vector<std::int64_t>>
            primeList(myRange, std::vector<std::int64_t>());

        double* tempEuler = nullptr;
        std::vector<std::int64_t> tempVec;
        MotleyMain(myMin, myMax, IsEuler, tempEuler,
                   tempVec, primeList, nThreads, maxThreads);

        cpp11::sexp myList = Rf_allocVector(VECSXP, myRange);

        for (std::size_t i = 0; i < myRange; ++i) {
            SET_VECTOR_ELT(myList, i, GetInt64Vec(primeList[i]));
        }

        if (keepNames) {
            ++numUnprotects;
            SetDblNames(myList, myRange,
                        static_cast<double>(myMin), myMax);
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

    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    const bool IsEuler = CleanConvert::convertFlag(RIsEuler, "IsEuler");
    const std::string namedObject = (IsEuler) ? "namedVector" : "namedList";
    bool IsNamed = CleanConvert::convertFlag(RNamed, namedObject);
    CleanConvert::convertPrimitive(Rb1, bound1, VecType::Numeric, "bound1");

    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, VecType::Numeric, "bound2");
    }

    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }

    if (myMax < 2) {
        if (IsEuler) {
            cpp11::sexp res = Rf_allocVector(INTSXP, 1);
            INTEGER(res)[0] = 1;

            if (IsNamed) {
                Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            }

            return res;
        } else {
            cpp11::sexp res = Rf_allocVector(VECSXP, 1);
            SET_VECTOR_ELT(res, 0, Rf_allocVector(INTSXP, 0));

            if (IsNamed) {
                Rf_setAttrib(res, R_NamesSymbol, Rf_mkString("1"));
            }

            return res;
        }
    }

    if (!Rf_isNull(RNumThreads)) {
        CleanConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    if (myMax > std::numeric_limits<int>::max()) {
        std::int64_t intMin = static_cast<std::int64_t>(myMin);
        return GlueDblMotley(intMin, myMax, IsEuler,
                             IsNamed, nThreads, maxThreads);
    } else {
        int intMin = static_cast<int>(myMin);
        int intMax = static_cast<int>(myMax);
        return GlueIntMotley(intMin, intMax, IsEuler,
                             IsNamed, nThreads, maxThreads);
    }
}
