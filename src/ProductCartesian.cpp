#include "cpp11/logicals.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"

#include "SetUpUtils.h"
#include "RMatrix.h"
#include <functional>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <thread>

[[cpp11::register]]
SEXP ExpandGridCpp(cpp11::list RList, SEXP Rlow, SEXP Rhigh,
                   SEXP RNumThreads, SEXP RmaxThreads) {

    const int nCols = Rf_length(RList);
    std::vector<int> IsFactor(nCols);
    std::vector<int> lenGrps(nCols);

    for (int i = 0; i < nCols; ++i) {
        if (Rf_isFactor(RList[i])) {
            IsFactor[i] = 1;
        } else {
            IsFactor[i] = 0;
        }

        lenGrps[i] = Rf_length(RList[i]);
    }

    const int sumLength = std::accumulate(
        lenGrps.begin(), lenGrps.end(), 0
    );

    std::vector<std::vector<int>> myVec(nCols);
    std::vector<int> typeCheck(N_TYPES, 0);

    cpp11::writable::strings charVec(sumLength);
    std::vector<Rcomplex> cmplxVec(sumLength);
    std::vector<Rbyte> rawVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    VecType myType = VecType::Integer;

    for (int i = 0, strt = 0; i < nCols; ++i) {
        switch(TYPEOF(RList[i])) {
            case INTSXP : {
                if (IsFactor[i]) {
                    typeCheck[tFac] = 1;
                } else {
                    typeCheck[tInt] = 1;
                }

                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), intVec.begin() + strt);
                myType = VecType::Integer;
                break;
            } case LGLSXP : {
                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), boolVec.begin() + strt);
                typeCheck[tLog] = 1;
                myType = VecType::Logical;
                break;
            } case CPLXSXP : {
                std::vector<Rcomplex> temp =
                    CppConvert::GetVec<Rcomplex>(RList[i]);
                std::copy(temp.begin(), temp.end(), cmplxVec.begin() + strt);
                typeCheck[tCpx] = 1;
                myType = VecType::Complex;
                break;
            } case RAWSXP : {
                std::vector<Rbyte> temp = CppConvert::GetVec<Rbyte>(RList[i]);
                std::copy(temp.begin(), temp.end(), rawVec.begin() + strt);
                typeCheck[tRaw] = 1;
                myType = VecType::Raw;
                break;
            } case REALSXP : {
                std::vector<double> temp =
                    CppConvert::GetVec<double>(RList[i]);
                std::copy(temp.begin(), temp.end(), dblVec.begin() + strt);
                typeCheck[tDbl] = 1;
                myType = VecType::Numeric;
                break;
            } case STRSXP : {
                for (int j = 0; j < lenGrps[i]; ++j) {
                    charVec[strt + j] = STRING_ELT(RList[i], j);
                }

                typeCheck[tStr] = 1;
                myType = VecType::Character;
                break;
            }
        }

        std::vector<int> idx(lenGrps[i]);
        std::iota(idx.begin(), idx.end(), strt);

        myVec[i] = idx;
        strt += lenGrps[i];
    }

    int mySum = std::accumulate(typeCheck.cbegin(), typeCheck.cend(), 0);

    // We need to check to see if there is overlap in factor levels
    if (typeCheck[tFac] && mySum == 1) {
        std::vector<std::string> testLevels;

        for (int i = 0; i < nCols; ++i) {
            if (IsFactor[i]) {
                cpp11::strings facVec(Rf_getAttrib(RList[i], R_LevelsSymbol));
                std::vector<std::string> strVec;

                for (auto f: facVec) {
                    const std::string temp(CHAR(f));
                    strVec.push_back(temp);
                }

                if (testLevels.size()) {
                    if (strVec != testLevels) {
                        ++mySum;
                        break;
                    }
                } else {
                    testLevels = strVec;
                }
            }
        }
    }

    bool IsDF = (mySum > 1) ? true : false;
    int nRows = 0;
    int nThreads = 1;
    int maxThreads = 1;

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    const double computedRows = CartesianCount(lenGrps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_class computedRowsMpz;

    if (IsGmp) {
        CartesianCountGmp(computedRowsMpz, lenGrps);
    }

    double lower = 0;
    double upper = 0;
    bool Parallel = false;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz, upperMpz, computedRowsMpz, computedRows);

    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                  lowerMpz, lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RNumThreads, limit);
    std::vector<int> startZ = nthProduct(lower, lenGrps);

    const int maxLen = *std::max_element(lenGrps.begin(), lenGrps.end());
    std::vector<int> mat_idx(maxLen * nCols);

    // transposing myVec so that each row represents
    // each element of the given list
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < lenGrps[i]; ++j) {
            mat_idx[i + j * nCols] = myVec[i][j];
        }
    }

    // Transform lenGrps to be used in nextProduct
    for (auto &v_i: lenGrps) {
        v_i = nCols * (v_i - 1);
    }

    cpp11::sexp res = GlueProdCart(
        mat_idx, typeCheck, IsFactor, RList, intVec, dblVec,
        boolVec, cmplxVec, rawVec, charVec, lenGrps, startZ,
        nRows, nCols, IsDF, nThreads, Parallel
    );

    return res;
}
