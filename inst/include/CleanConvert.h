#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

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
    Integer = 1,
    Numeric = 2,
    Logical = 3,
    Character = 4,
    Complex = 5,
    Raw = 6
};

// PartitionEsque = 2: Can't be reduced to an integer partition but still has similarities
// to the more general subset sum problem. E.g. v = rnorm(20, mean = 10.5), m = 4,
// rep = TRUE, tar = c(11.005, 11.15), comparisonFun = c(">", "<"), constraintFun = "mean"

// PartGeneral = 3: Occurs when non-standard input can be reduced to a general integer
// partition: E.g. v = seq(200, 300, 5), tar = 1200, m = 4, rep = TRUE

enum class PartitionType {
    NotPartition = 1,
    PartitonEsque = 2,
    PartGeneral = 3,
    PartTraditional = 4, // Get all partitions. E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    PartTradNoZero = 5, // E.g. tar = 20 startZ = c(1, 1, 1, 1, 15)
    PartDstctStdAll = 6, // Get all distinct partitions (0 can repeat) E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    PartDstctShort = 7, // Case where startZ doesn't maximize width. E.g. tar = 20 startZ = c(0, 0, 20)
    PartDstctSpecial = 8, // Case where startZ doesn't maximize 0's. E.g. tar = 20 startZ = c(0, 0, 1, 2, 17)
    PartDstctOneZero = 9, // Similar to above but can occur when IsMult = FALSE. E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
    PartDstctNoZero = 10, // E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
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
