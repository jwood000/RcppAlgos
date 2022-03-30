#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

#include <limits>
#include <vector>
#include <string>
#include <gmp.h>

// std::pow(2, 53) - 1 is not constexpr as of C++11
constexpr double Significand53 = 9007199254740991.0;

// Based off the internal limitations of sample, we
// cannot utilize the full range of 53-bit significand
// precision. Here is the condition from `do_dample2`:
//     if (!R_FINITE(dn) || dn < 0 || dn > 4.5e15 || (k > 0 && dn == 0))
// Here is the source (line ~2100):
//     https://github.com/wch/r-source/blob/trunk/src/main/unique.c
constexpr double SampleLimit = 4500000000000000.0;

enum class VecType {
    Integer   = 1,
    Numeric   = 2,
    Logical   = 3,
    Character = 4,
    Complex   = 5,
    Raw       = 6
};

namespace CleanConvert {

    bool CheckNA(double val, VecType myType);

    template <typename T>
    std::vector<T> GetNumVec(SEXP Rv);

    SEXP GetCount(bool IsGmp, const mpz_t computedRowsMpz,
                  double computedRows);

    bool convertFlag(SEXP boolInput, const std::string &nameOfBool);

    template <typename T>
    void convertPrimitive(SEXP input, T &result, VecType myType,
                          const std::string &nameOfObject,
                          bool numOnly = true, bool checkWhole = true,
                          bool negPoss = false, bool decimalFraction = false);

    template <typename T>
    void convertVector(SEXP input, std::vector<T> &result,
                       VecType myType, const std::string &nameOfObject,
                       bool numOnly = true, bool checkWhole = true,
                       bool negPoss = false);
}

#endif
