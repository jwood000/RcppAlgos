#include "NumbersUtils/Eratosthenes.h"
#include "CleanConvert.h"

[[cpp11::register]]
SEXP PrimeSieveCpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads,
                   SEXP RmaxCores, SEXP RmaxThreads) {

    double bound1;
    double bound2;

    std::int_fast64_t myMin;
    std::int_fast64_t myMax;

    int maxCores = 1;
    int nThreads = 1;
    int maxThreads = 1;

    CleanConvert::convertPrimitive(RmaxCores, maxCores,
                                   VecType::Integer, "maxCores");
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    CleanConvert::convertPrimitive(Rb1, bound1, VecType::Numeric,
                                   "bound1", true, false);

    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, VecType::Numeric,
                                       "bound2", true, false);
    }

    if (bound1 > bound2) {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound1));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound2));
    } else {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound2));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound1));
    }

    if (myMax <= 1) {
        return Rf_allocVector(INTSXP, 0);
    }

    if (myMin <= 2) myMin = 1;

    if (myMin == myMax) {
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
        CleanConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    int numPrimes = 0;
    std::vector<int> runningCount;

    runningCount.push_back(0);
    int numSects = nThreads;
    bool Parallel = false;

    if (myMax > std::numeric_limits<int>::max()) {
        std::vector<std::vector<double>> primeList(numSects,
                                                   std::vector<double>());
        std::vector<double> tempPrimes;

        PrimeSieve::PrimeSieveMain(primeList, tempPrimes, myMin, myMax,
                                   Parallel, nThreads, maxThreads, maxCores);

        if (Parallel) {
            for (int i = 0; i < numSects; ++i) {
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
        std::vector<std::vector<int>> primeList(numSects,
                                                std::vector<int>());
        std::vector<int> tempPrimes;

        PrimeSieve::PrimeSieveMain(primeList, tempPrimes, myMin, myMax,
                                   Parallel, nThreads, maxThreads, maxCores);

        if (Parallel) {
            for (int i = 0; i < numSects; ++i) {
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
