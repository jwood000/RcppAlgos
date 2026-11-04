#include "NumbersUtils/PrimeCount.h"
#include "CppConvert.h"
#include <limits>

[[cpp11::register]]
SEXP PrimeCountCpp(SEXP Rn, SEXP RNumThreads, SEXP RmaxThreads) {
    double dblNum;
    CppConvert::convertPrimitive(Rn, dblNum, VecType::Numeric, "n");
    const std::int64_t n = static_cast<std::int64_t>(dblNum);

    int nThreads = 1;
    int maxThreads = 1;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    if (n < 100000) {
        if (n < 10) {
            if (n == 1)
                return Rf_ScalarInteger(0);
            else if (n == 2)
                return Rf_ScalarInteger(1);
            else if (n < 5)
                return Rf_ScalarInteger(2);
            else if (n < 7)
                return Rf_ScalarInteger(3);
            else
                return Rf_ScalarInteger(4);
        }

        return Rf_ScalarInteger(static_cast<int>(PrimeCounting::PiPrime(n)));
    }

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    std::int64_t result = PrimeCounting::PiLegendre(n, nThreads, maxThreads);

    if (result > std::numeric_limits<int>::max()) {
        return Rf_ScalarReal(static_cast<double>(result));
    } else {
        return Rf_ScalarInteger(static_cast<int>(result));
    }
}
