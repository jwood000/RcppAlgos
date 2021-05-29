#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

#include <vector>
#include <string>

constexpr double Significand53 = 9007199254740991.0;

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

    bool convertLogical(SEXP boolInput, const std::string &nameOfBool);

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
