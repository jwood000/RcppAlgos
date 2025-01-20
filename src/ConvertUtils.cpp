#include "CppConvert/ConvertUtils.h"
#include <numeric>

namespace CppConvert {

    template <typename T>
    void SetNames(SEXP res, T myMin, T myMax) {
        cpp11::writable::r_vector<T> myNames((myMax - myMin) + 1);
        std::iota(myNames.begin(), myNames.end(), myMin);
        Rf_setAttrib(res, R_NamesSymbol, myNames);
    }

    template <typename T>
    void SetNames(SEXP res, const std::vector<T> &myNums) {
        cpp11::writable::r_vector<T> myNames(myNums);
        Rf_setAttrib(res, R_NamesSymbol, myNames);
    }

    bool CheckNA(double val, VecType myType) {
        if (myType == VecType::Integer) {
            return (ISNAN(val) || val == NA_INTEGER);
        } else {
            return ISNAN(val);
        }
    }

    template <>
    std::vector<Rcomplex> GetVec(SEXP Rv) {
        Rcomplex* cmplxRv = COMPLEX(Rv);
        std::vector<Rcomplex> v(cmplxRv, cmplxRv + Rf_length(Rv));
        return v;
    }

    template <>
    std::vector<Rbyte> GetVec(SEXP Rv) {
        Rbyte* rawRv = RAW(Rv);
        std::vector<Rbyte> v(rawRv, rawRv + Rf_length(Rv));
        return v;
    }

    template <typename T>
    std::vector<T> GetVec(SEXP Rv) {
        std::vector<T> v;
        const int len = Rf_length(Rv);

        if (len) {
            switch(TYPEOF(Rv)) {
                case LGLSXP : {
                    int* boolRv = LOGICAL(Rv);
                    v.assign(boolRv, boolRv + len);
                    break;
                } case REALSXP : {
                    double* dblRv = REAL(Rv);
                    v.assign(dblRv, dblRv + len);
                    break;
                } case INTSXP : {
                    int* intRv = INTEGER(Rv);
                    v.assign(intRv, intRv + len);
                    break;
                }
            }
        }

        return v;
    }

    int rawExport(char* raw, mpz_class value, std::size_t totals) {

        std::memset(raw, 0, totals);

        int* r = (int*)raw;
        r[0] = totals / intSize - 2;

        r[1] = static_cast<int>(mpz_sgn(value.get_mpz_t()));
        mpz_export(&r[2], 0, 1, intSize, 0, 0, value.get_mpz_t());

        return totals;
    }

    void QuickSort(std::vector<mpz_class> &arr, int left,
                   int right, std::vector<std::size_t> &lens) {

        int i = left;
        int j = right;
        mpz_class pivot;

        int mid = (left + right) / 2;
        pivot = arr[mid];

        /* partition */
        while (i <= j) {
            while (cmp(arr[i], pivot) < 0)
                ++i;

            while (j >= 0 && cmp(arr[j], pivot) > 0)
                --j;

            if (i <= j && j >= 0) {
                swap(arr[i], arr[j]);
                std::swap(lens[i], lens[j]);
                ++i;
                --j;
            }
        }

        /* recursion */
        if (left < j) QuickSort(arr, left, j, lens);
        if (i < right) QuickSort(arr, i, right, lens);
    }

    SEXP GetCount(bool IsGmp, mpz_class numMpz, double numDbl) {

        if (IsGmp) {
            const std::size_t sizeNum = intSize *
                (2 + (mpz_sizeinbase(numMpz.get_mpz_t(), 2) + numb - 1) / numb);
            const std::size_t size = intSize + sizeNum;

            cpp11::sexp ans = Rf_allocVector(RAWSXP, size);
            char* rPos = (char*) RAW(ans);
            ((int*) rPos)[0] = 1; // first int is vector-size-header

            // current position in rPos[] (starting after vector-size-header)
            rawExport(&rPos[intSize], numMpz, sizeNum);
            Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
            return(ans);
        } else {
            if (numDbl > std::numeric_limits<int>::max()) {
                return Rf_ScalarReal(numDbl);
            } else {
                return Rf_ScalarInteger(static_cast<int>(numDbl));
            }
        }
    }
}

template void CppConvert::SetNames(SEXP, double, double);
template void CppConvert::SetNames(SEXP, int, int);

template void CppConvert::SetNames(SEXP, const std::vector<double>&);
template void CppConvert::SetNames(SEXP, const std::vector<int>&);

template std::vector<double> CppConvert::GetVec(SEXP);
template std::vector<int> CppConvert::GetVec(SEXP);
