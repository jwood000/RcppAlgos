#include "NumbersUtils/Eratosthenes.h"
#include "NumbersUtils/PrimeCount.h"
#include "CppConvert.h"
#include "cpp11/as.hpp"

[[cpp11::register]]
SEXP PrimeSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads,
                   SEXP RmaxCores, SEXP RmaxThreads, SEXP RCountOnly) {

    const bool countOnly = CppConvert::convertFlag(RCountOnly, "countOnly");

    double bound1;
    double bound2;

    std::int_fast64_t myMin;
    std::int_fast64_t myMax;

    int maxCores = 1;
    int nThreads = 1;
    int maxThreads = 1;

    CppConvert::convertPrimitive(RmaxCores, maxCores,
                                   VecType::Integer, "maxCores");
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    CppConvert::convertPrimitive(Rb1, bound1, VecType::Numeric,
                                   "bound1", true, false);

    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CppConvert::convertPrimitive(Rb2, bound2, VecType::Numeric,
                                       "bound2", true, false);
    }

    if (bound1 > bound2) {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound1));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound2));
    } else {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound2));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound1));
    }

    if (!countOnly && myMax <= 1) {
        return Rf_allocVector(INTSXP, 0);
    }

    if (myMin <= 2) myMin = 1;

    if (!countOnly && myMin == myMax) {
        if (myMax % 2) {
            ++myMax;
        } else {
            if (myMax > std::numeric_limits<int>::max()) {
                return Rf_allocVector(REALSXP, 0);
            } else {
                return Rf_allocVector(INTSXP, 0);
            }
        }
    }

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    if (countOnly) {
        // PiPrimeRange handles empty and singleton intervals directly. It also
        // normalizes the budgets before calling the shared sieve launcher.
        const std::int64_t count = PrimeCounting::PiPrimeRange(
            myMin, myMax, nThreads, maxThreads, maxCores
        );

        if (count > std::numeric_limits<int>::max()) {
            return cpp11::as_sexp(static_cast<double>(count));
        } else {
            return cpp11::as_sexp(static_cast<int>(count));
        }
    }

    int numPrimes = 0;
    std::vector<int> runningCount;

    runningCount.push_back(0);
    int numSects = nThreads;
    bool Parallel = false;

    if (myMax > std::numeric_limits<int>::max()) {
        std::vector<std::vector<double>> primeList;
        std::vector<double> tempPrimes;

        PrimeSieve::PrimeSieveMain(primeList, tempPrimes, myMin, myMax,
                                   Parallel, nThreads, maxThreads, maxCores);

        if (Parallel) {
            const std::size_t numSects = primeList.size();

            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }

            cpp11::sexp res = Rf_allocVector(REALSXP, numPrimes);
            double* primes = REAL(res);

            for (int i = 0; i < numSects; ++i) {
                std::copy(primeList[i].cbegin(), primeList[i].cend(),
                          primes + runningCount[i]);
            }

            return res;
        } else {
            cpp11::sexp primes = Rf_allocVector(REALSXP, tempPrimes.size());
            std::copy(tempPrimes.cbegin(), tempPrimes.cend(), REAL(primes));
            return primes;
        }
    } else {
        std::vector<std::vector<int>> primeList;
        std::vector<int> tempPrimes;

        PrimeSieve::PrimeSieveMain(primeList, tempPrimes, myMin, myMax,
                                   Parallel, nThreads, maxThreads, maxCores);

        if (Parallel) {
            const std::size_t numSects = primeList.size();

            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }

            cpp11::sexp res = Rf_allocVector(INTSXP, numPrimes);
            int* primes = INTEGER(res);

            for (int i = 0; i < numSects; ++i) {
                std::copy(primeList[i].cbegin(), primeList[i].cend(),
                          primes + runningCount[i]);
            }

            return res;
        } else {
            cpp11::sexp primes = Rf_allocVector(INTSXP, tempPrimes.size());
            std::move(tempPrimes.cbegin(), tempPrimes.cend(), INTEGER(primes));
            return primes;
        }
    }
}
