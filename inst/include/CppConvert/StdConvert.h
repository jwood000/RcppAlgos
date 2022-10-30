#pragma once

#include "cpp11/protect.hpp"
#include "ConvertUtils.h"

namespace CppConvert {

    bool convertFlag(SEXP boolInput, const std::string &nameOfBool);

    template <typename T>
    void convertPrimitive(SEXP input, T &result, VecType myType,
                          const std::string &nameOfObject, bool numOnly = true,
                          bool checkWhole = true, bool negPoss = false,
                          bool decimalFraction = false);

    template <typename T>
    void convertVector(SEXP input, std::vector<T> &result, VecType myType,
                       const std::string &nameOfObject, bool numOnly = true,
                       bool checkWhole = true, bool negPoss = false);
}
