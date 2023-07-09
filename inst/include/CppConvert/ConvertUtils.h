#pragma once

#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/sexp.hpp"

#include "Constants.h"
#include <gmpxx.h>
#include <vector>
#include <string>
#include <cmath>

namespace CppConvert {

    template <typename T>
    void SetNames(SEXP res, T myMin, T myMax);

    template <typename T>
    void SetNames(SEXP res, const std::vector<T> &myNums);

    bool CheckNA(double val, VecType myType);

    template <typename T>
    std::vector<T> GetNumVec(SEXP Rv);

    int rawExport(char* raw, mpz_class value, std::size_t totals);
    void QuickSort(std::vector<mpz_class> &arr, int left,
                   int right, std::vector<std::size_t> &lens);
    SEXP GetCount(bool IsGmp, mpz_class numMpz, double numDbl);
}
