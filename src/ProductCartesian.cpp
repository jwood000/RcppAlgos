#include "cpp11/logicals.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"

#include "SetUpUtils.h"
#include <algorithm>
#include <cstdint>
#include <numeric>

double CartesianCount(const std::vector<int> &lenGrps) {
    return std::accumulate(lenGrps.begin(), lenGrps.end(),
                           1.0, std::multiplies<double>());
}

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps) {

    result = 1;

    for (auto len: lenGrps) {
        result *= len;
    }
}

void GetCharOutput(cpp11::writable::strings_matrix<> &mat,
                   const std::vector<std::vector<int>> &idx,
                   const cpp11::strings &charVec,
                   int nCols, int nRows) {

    const int lastCol = nCols - 1;

    for (int i = 0, preSecLen = nRows; i < lastCol; ++i) {
        for (int j = 0, q = 0, grpSize = idx[i].size(),
             secLen = preSecLen / grpSize; j < nRows; j += secLen) {

            // Benchmarks show that this set up is about 30%
            // faster than using cpp11::sexp
            SEXP val = PROTECT(STRING_ELT(charVec, idx[i][q]));

            for (int k = 0, m_idx = i * nRows + j; k < secLen; ++k) {
                SET_STRING_ELT(mat, m_idx + k, val);
            }

            UNPROTECT(1);
            ++q;
            q %= grpSize;
        }

        preSecLen /= idx[i].size();
    }

    const int lastGrpSize = idx.back().size();

    for (int j = 0, q = charVec.size() - lastGrpSize,
         secLen = lastGrpSize; j < nRows; j += secLen) {
        for (int k = 0, m_idx = lastCol * nRows + j; k < secLen; ++k) {
            SET_STRING_ELT(mat, m_idx + k, STRING_ELT(charVec, q + k));
        }
    }
}

template <typename matType, typename T>
void GetPureOutput(matType &mat,
                   const std::vector<std::vector<int>> &idx,
                   const std::vector<T> &standardVec,
                   int nCols, int nRows) {

    const int lastCol = nCols - 1;

    for (int i = 0, preSecLen = nRows; i < lastCol; ++i) {
        for (int j = 0, q = 0, grpSize = idx[i].size(),
             secLen = preSecLen / grpSize; j < nRows; j += secLen) {

            auto&& val = standardVec[idx[i][q]];

            for (int k = 0; k < secLen; ++k) {
                mat(j + k, i) = val;
            }

            ++q;
            q %= grpSize;
        }

        preSecLen /= idx[i].size();
    }

    const int lastGrpSize = idx.back().size();
    std::vector<T> lastVec(lastGrpSize);
    std::copy(standardVec.end() - lastGrpSize,
              standardVec.end(), lastVec.begin());

    for (int j = 0, secLen = lastGrpSize; j < nRows; j += secLen) {
        for (int k = 0; k < secLen; ++k) {
            mat(j + k, lastCol) = lastVec[k];
        }
    }
}


template <typename cpp11Type, typename T>
void PoulateColumn(cpp11Type &res,
                   const std::vector<int> &idx,
                   const std::vector<T> &poolVec,
                   int nCols, int nRows, int preSecLen, int i) {

    if (i < (nCols - 1)) {
        for (int j = 0, q = 0, grpSize = idx.size(),
             secLen = preSecLen / grpSize; j < nRows;) {

            auto&& val = poolVec[idx[q]];

            for (int k = 0; k < secLen; ++k, ++j) {
                res[j] = val;
            }

            ++q;
            q %= grpSize;
        }
    } else {
        const int lastGrpSize = idx.size();
        std::vector<T> lastVec(lastGrpSize);
        std::copy(poolVec.end() - lastGrpSize,
                  poolVec.end(), lastVec.begin());

        for (int j = 0, secLen = lastGrpSize; j < nRows;) {
            for (int k = 0; k < secLen; ++k, ++j) {
                res[j] = lastVec[k];
            }
        }
    }
}

SEXP GlueProdCart(
    const std::vector<std::vector<int>> &idx,
    const std::vector<int> &typeCheck, const std::vector<int> &IsFactor,
    const cpp11::list &RList, const std::vector<int> &intVec,
    const std::vector<double> &dblVec, const std::vector<int> &boolVec,
    const cpp11::strings &charVec, int nRows, int nCols, bool IsDF
) {

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);

        for (int i = 0, preSecLen = nRows; i < nCols; ++i) {
            switch (TYPEOF(RList[i])) {
                case INTSXP: {
                    cpp11::writable::integers intSexpVec(nRows);

                    if (IsFactor[i]) {
                        PoulateColumn(intSexpVec, idx[i], intVec,
                                      nCols, nRows, preSecLen, i);
                        SetFactorClass(intSexpVec, RList[i]);
                    } else {
                        PoulateColumn(intSexpVec, idx[i], intVec,
                                      nCols, nRows, preSecLen, i);
                    }

                    DataFrame[i] = intSexpVec;
                    break;
                } case LGLSXP: {
                    cpp11::writable::logicals boolSexpVec(nRows);
                    PoulateColumn(boolSexpVec, idx[i], boolVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = boolSexpVec;
                    break;
                } case REALSXP: {
                    cpp11::writable::doubles dblSexpVec(nRows);
                    PoulateColumn(dblSexpVec, idx[i], dblVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = dblSexpVec;
                    break;
                } case STRSXP: {
                    cpp11::writable::strings sexpVec(nRows);

                    if (i < (nCols - 1)) {
                        for (int j = 0, q = 0, grpSize = idx[i].size(),
                             secLen = preSecLen / grpSize; j < nRows;) {

                            // Benchmarks show that this set up is about 30%
                            // faster than using cpp11::sexp
                            SEXP val = PROTECT(STRING_ELT(charVec, idx[i][q]));

                            for (int k = 0; k < secLen; ++k, ++j) {
                                SET_STRING_ELT(sexpVec, j, val);
                            }

                            UNPROTECT(1);
                            ++q;
                            q %= grpSize;
                        }
                    } else {
                        const int lastGrpSize = idx.back().size();

                        for (int j = 0, q = charVec.size() - lastGrpSize,
                             secLen = lastGrpSize; j < nRows;) {
                            for (int k = 0; k < secLen; ++k, ++j) {
                                SET_STRING_ELT(
                                    sexpVec, j, STRING_ELT(charVec, q + k)
                                );
                            }
                        }
                    }

                    DataFrame[i] = sexpVec;
                    break;
                }
            }

            preSecLen /= idx[i].size();
        }

        DataFrame.attr("row.names") = {NA_INTEGER, -nRows};
        DataFrame.names() = RList.names();
        DataFrame.attr("class") = "data.frame";
        return DataFrame;
    } else {
        if (typeCheck[tInt]) {
            cpp11::writable::integers_matrix<> intMat(nRows, nCols);
            GetPureOutput(intMat, idx, intVec, nCols, nRows);
            return intMat;
        } else if (typeCheck[tFac]) {
            cpp11::writable::integers_matrix<> facMat(nRows, nCols);
            GetPureOutput(facMat, idx, intVec, nCols, nRows);
            SetFactorClass(facMat, RList[0]);
            return facMat;
        } else if (typeCheck[tLog]) {
            cpp11::writable::logicals_matrix<> boolMat(nRows, nCols);
            GetPureOutput(boolMat, idx, boolVec, nCols, nRows);
            return boolMat;
        } else if (typeCheck[tDbl]) {
            cpp11::writable::doubles_matrix<> dblMat(nRows, nCols);
            GetPureOutput(dblMat, idx, dblVec, nCols, nRows);
            return dblMat;
        } else {
            cpp11::writable::strings_matrix<> charMat(nRows, nCols);
            GetCharOutput(charMat, idx, charVec, nCols, nRows);
            return charMat;
        }
    }
}

[[cpp11::register]]
SEXP ExpandGridCpp(cpp11::list RList) {

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
    std::vector<int> typeCheck(5, 0);

    cpp11::writable::strings charVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    for (int i = 0, strt = 0; i < nCols; ++i) {
        switch(TYPEOF(RList[i])) {
            case INTSXP: {
                if (IsFactor[i]) {
                    typeCheck[tFac] = 1;
                } else {
                    typeCheck[tInt] = 1;
                }

                std::vector<int> temp = CppConvert::GetNumVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), intVec.begin() + strt);
                break;
            } case LGLSXP: {
                std::vector<int> temp = CppConvert::GetNumVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), boolVec.begin() + strt);
                typeCheck[tLog] = 1;
                break;
            } case REALSXP: {
                std::vector<double> temp =
                    CppConvert::GetNumVec<double>(RList[i]);
                std::copy(temp.begin(), temp.end(), dblVec.begin() + strt);
                typeCheck[tDbl] = 1;
                break;
            } case STRSXP: {
                for (int j = 0; j < lenGrps[i]; ++j) {
                    charVec[strt + j] = STRING_ELT(RList[i], j);
                }

                typeCheck[tStr] = 1;
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
    int nRows = CartesianCount(lenGrps);

    cpp11::sexp res = GlueProdCart(
        myVec, typeCheck, IsFactor, RList, intVec,
        dblVec, boolVec, charVec, nRows, nCols, IsDF
    );

    return res;
}