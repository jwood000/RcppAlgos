#include "NumbersUtils/Eratosthenes.h"
#include "CartesianContainer.h"
#include "ComboCartesian.h"
#include "CleanConvert.h"
#include "SetUpUtils.h"
#include <unordered_map>
#include <numeric>  // std::accumulate

enum rcppType {
    tInt = 0,
    tDbl = 1,
    tStr = 2,
    tLog = 3,
    tFac = 4
};

void convertToString(std::vector<std::string> &tempVec,
                     SEXP ListElement, rcppType &typePass, bool bFac) {

    const int len = Rf_length(ListElement);

    switch(TYPEOF(ListElement)) {
        case INTSXP: {
            if (bFac) {
                SEXP facVec = PROTECT(Rf_coerceVector(ListElement, INTSXP));
                SEXP myLevels = PROTECT(Rf_getAttrib(facVec, R_LevelsSymbol));
                std::vector<std::string> strVec;

                for (int i = 0, faclen = Rf_length(myLevels); i < faclen; ++i) {
                    strVec.push_back(CHAR(STRING_ELT(myLevels, i)));
                }

                UNPROTECT(2);
                typePass = tFac;

                for (auto v: strVec) {
                    bool isNum = !v.empty() &&
                                 v.find_first_not_of("0123456789.") ==
                                 std::string::npos;

                    if (isNum) {
                        double n = std::atof(v.c_str());
                        if (std::floor(n) != n) v += "dbl";
                    } else {
                        v += "str";
                    }

                    tempVec.push_back(v);
                }
            } else {
                int* intPtr = INTEGER(ListElement);
                std::vector<int> intVec(intPtr, intPtr + len);
                typePass = tInt;

                for (auto v: intVec) {
                    tempVec.push_back(std::to_string(v));
                }
            }

            break;
        } case LGLSXP: {
            int* intPtr = INTEGER(ListElement);
            std::vector<int> intVec(intPtr, intPtr + len);
            typePass = tLog;

            for (auto v: intVec) {
                std::string r = std::to_string(v);
                tempVec.push_back(r + "log");
            }

            break;
        } case REALSXP: {
            double* dblPtr = REAL(ListElement);
            std::vector<double> dblVec(dblPtr, dblPtr + len);
            typePass = tDbl;

            for (auto v: dblVec) {
                std::string r = std::to_string(v);

                while (r.back() == '0') {
                    r.pop_back();
                }

                if (std::floor(v) != v) r += "dbl";
                tempVec.push_back(r);
            }

            break;
        } case STRSXP: {
            std::vector<std::string> strVec;
            typePass = tStr;

            for (int i = 0; i < len; ++i) {
                const std::string temp(CHAR(STRING_ELT(ListElement, i)));
                strVec.push_back(temp);
            }

            for (auto v: strVec) {
                bool isNum = !v.empty() &&
                    v.find_first_not_of("0123456789.") ==
                    std::string::npos;

                if (isNum) {
                    double n = std::atof(v.c_str());
                    if (std::floor(n) != n) v += "dbl";
                } else {
                    v += "str";
                }

                tempVec.push_back(v);
            }

            break;
        }
    }
}

void AddNames(SEXP res, SEXP RList) {
    if (!Rf_isNull(Rf_getAttrib(RList, R_NamesSymbol))) {
        SEXP myNames = PROTECT(Rf_getAttrib(RList, R_NamesSymbol));
        SEXP dimNames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimNames, 1, myNames);
        Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
        UNPROTECT(2);
    }
}

void GetCharOutput(SEXP res, SEXP RList,
                   const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   SEXP charVec, int nCols, int nRows) {

    for (int i = 0, n1 = nCols - 1, row = 0, m_idx = 0,
         baseSize = lenGrps.size(); i < baseSize; ++i, m_idx = row) {

        for (int j = 0, c_idx = i * n1, grpSize = lenGrps[i];
             j < n1; ++j, ++c_idx, m_idx += nRows) {

            SEXP comb = PROTECT(STRING_ELT(charVec, cartCombs[c_idx]));

            for (int k = 0; k < grpSize; ++k) {
                SET_STRING_ELT(res, m_idx + k, comb);
            }

            UNPROTECT(1);
        }

        for (int k = 0; k < lenGrps[i]; ++k, ++row) {
            SET_STRING_ELT(res, m_idx + k,
                           STRING_ELT(charVec, lastCol[row]));
        }
    }

    AddNames(res, RList);
}

template <typename T>
void GetPureOutput(T* mat, SEXP res, SEXP RList,
                   const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const T* standardVec, int nCols, int nRows) {

    for (int i = 0, n1 = nCols - 1, row = 0, m_idx = 0,
         baseSize = lenGrps.size(); i < baseSize; ++i, m_idx = row) {

        for (int j = 0, c_idx = i * n1, grpSize = lenGrps[i];
             j < n1; ++j, ++c_idx, m_idx += nRows) {

            auto&& comb = standardVec[cartCombs[c_idx]];

            for (int k = 0; k < grpSize; ++k) {
                mat[m_idx + k] = comb;
            }
        }

        for (int k = 0; k < lenGrps[i]; ++k, ++row) {
            mat[m_idx + k] = standardVec[lastCol[row]];
        }
    }

    AddNames(res, RList);
}

SEXP GlueComboCart(const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const std::vector<std::vector<int>> &facList,
                   const std::vector<int> &typeCheck,
                   const std::vector<int> &IsFactor, SEXP RList,
                   int* intVec, int* boolVec, double* dblVec,
                   SEXP charVec, int nRows, int nCols, bool IsDF) {

    if (IsDF) {
        SEXP DataFrame = PROTECT(Rf_allocVector(VECSXP, nCols));
        int numProtects = 1;

        for (int i = 0, facInd = 0; i < nCols; ++i) {
            switch (TYPEOF(VECTOR_ELT(RList, i))) {
                case INTSXP: {
                    SEXP sexpVec = PROTECT(Rf_allocVector(INTSXP, nRows));
                    int* intSexpVec = INTEGER(sexpVec);
                    ++numProtects;

                    if (IsFactor[i]) {
                        const int size = facList[facInd].size();
                        SEXP facVec = PROTECT(Rf_allocVector(INTSXP, size));
                        int* intFacVec = INTEGER(facVec);
                        ++numProtects;

                        for (int j = 0; j < size; ++j) {
                            intFacVec[j] = facList[facInd][j];
                        }

                        if (i < (nCols - 1)) {
                            for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
                                 baseSize = lenGrps.size(); idx < baseSize;
                                 ++idx, row += n1) {

                                auto&& comb = intFacVec[cartCombs[row]];

                                for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                                    intSexpVec[j] = comb;
                                }
                            }
                        } else {
                            for (int j = 0; j < nRows; ++j) {
                                intSexpVec[j] = intFacVec[lastCol[j]];
                            }
                        }

                        SetFactorClass(sexpVec, VECTOR_ELT(RList, i));
                        ++facInd;
                    } else {
                        if (i < (nCols - 1)) {
                            for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
                                 baseSize = lenGrps.size(); idx < baseSize;
                                 ++idx, row += n1) {

                                auto&& comb = intVec[cartCombs[row]];

                                for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                                    intSexpVec[j] = comb;
                                }
                            }
                        } else {
                            for (int j = 0; j < nRows; ++j) {
                                intSexpVec[j] = intVec[lastCol[j]];
                            }
                        }
                    }

                    SET_VECTOR_ELT(DataFrame, i, sexpVec);
                    break;
                } case LGLSXP: {
                    SEXP sexpVec = PROTECT(Rf_allocVector(LGLSXP, nRows));
                    int* boolSexpVec = INTEGER(sexpVec);
                    ++numProtects;

                    if (i < (nCols - 1)) {
                        for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
                             baseSize = lenGrps.size(); idx < baseSize;
                             ++idx, row += n1) {

                            auto&& comb = boolVec[cartCombs[row]];

                            for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                                boolSexpVec[j] = comb;
                            }
                        }
                    } else {
                        for (int j = 0; j < nRows; ++j) {
                            boolSexpVec[j] = boolVec[lastCol[j]];
                        }
                    }

                    SET_VECTOR_ELT(DataFrame, i, sexpVec);
                    break;
                } case REALSXP: {
                    SEXP sexpVec = PROTECT(Rf_allocVector(REALSXP, nRows));
                    double* dblSexpVec = REAL(sexpVec);
                    ++numProtects;

                    if (i < (nCols - 1)) {
                        for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
                             baseSize = lenGrps.size(); idx < baseSize;
                             ++idx, row += n1) {

                            auto&& comb = dblVec[cartCombs[row]];

                            for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                                dblSexpVec[j] = comb;
                            }
                        }
                    } else {
                        for (int j = 0; j < nRows; ++j) {
                            dblSexpVec[j] = dblVec[lastCol[j]];
                        }
                    }

                    SET_VECTOR_ELT(DataFrame, i, sexpVec);
                    break;
                } case STRSXP: {
                    SEXP sexpVec = PROTECT(Rf_allocVector(STRSXP, nRows));
                    ++numProtects;

                    if (i < (nCols - 1)) {
                        for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
                             baseSize = lenGrps.size(); j < baseSize;
                             ++idx, row += n1) {

                            SEXP comb = PROTECT(STRING_ELT(charVec, cartCombs[row]));

                            for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                                SET_STRING_ELT(sexpVec, j, comb);
                            }

                            UNPROTECT(1);
                        }
                    } else {
                        for (int j = 0; j < nRows; ++j) {
                            SET_STRING_ELT(sexpVec, j,
                                           STRING_ELT(charVec, lastCol[j]));
                        }
                    }

                    SET_VECTOR_ELT(DataFrame, i, sexpVec);
                    break;
                }
            }
        }

        Rf_setAttrib(DataFrame, R_ClassSymbol, Rf_mkString("data.frame"));
        Rf_setAttrib(DataFrame, R_NamesSymbol,
                     Rf_getAttrib(RList, R_NamesSymbol));
        SEXP rownames = PROTECT(Rf_allocVector(INTSXP, nRows));
        int* intRows  = INTEGER(rownames);

        for (int i = 1; i <= nRows; ++i) {
            intRows[i - 1] = i;
        }

        Rf_setAttrib(DataFrame, R_RowNamesSymbol, rownames);
        UNPROTECT(numProtects + 1);
        return DataFrame;
    } else {
        if (typeCheck[tInt]) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, nCols));
            int* intMat = INTEGER(res);
            GetPureOutput(intMat, res, RList, cartCombs, lastCol,
                          lenGrps, intVec, nCols, nRows);
            UNPROTECT(1);
            return res;
        } else if (typeCheck[tFac]) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, nCols));
            int* intMat = INTEGER(res);
            GetPureOutput(intMat, res, RList, cartCombs, lastCol,
                          lenGrps, intVec, nCols, nRows);
            SetFactorClass(res, VECTOR_ELT(RList, 0));
            UNPROTECT(1);
            return res;
        } else if (typeCheck[tLog]) {
            SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, nRows, nCols));
            int* intMat = INTEGER(res);
            GetPureOutput(intMat, res, RList, cartCombs, lastCol,
                          lenGrps, boolVec, nCols, nRows);
            UNPROTECT(1);
            return res;
        } else if (typeCheck[tDbl]) {
            SEXP res = Rf_allocMatrix(REALSXP, nRows, nCols);
            double* dblMat = REAL(res);
            GetPureOutput(dblMat, res, RList, cartCombs, lastCol,
                          lenGrps, dblVec, nCols, nRows);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(STRSXP, nRows, nCols));
            GetCharOutput(res, RList, cartCombs, lastCol,
                          lenGrps, charVec, nCols, nRows);
            UNPROTECT(1);
            return res;
        }
    }
}

void getAtLeastNPrimes(std::vector<int> &primes,
                       std::size_t sumLength) {

    double limit = 100;
    std::size_t guess = PrimeSieve::EstimatePiPrime(1.0, limit);

    while (guess < (1.1 * sumLength)) {
        limit *= 2;
        guess = PrimeSieve::EstimatePiPrime(1.0, limit);
    }

    bool tempPar = false;
    int_fast64_t intMin = static_cast<int_fast64_t>(1);
    int_fast64_t intMax = static_cast<int_fast64_t>(limit);
    std::vector<std::vector<int>> tempList;

    PrimeSieve::PrimeSieveMain(tempList, primes, intMin, intMax, tempPar);
}

SEXP ComboGridCpp(SEXP RList, SEXP RIsRep) {

    int sumLength = 0;
    const int nCols = Rf_length(RList);
    std::vector<int> IsFactor(nCols);

    for (int i = 0; i < nCols; ++i) {
        if (Rf_isFactor(VECTOR_ELT(RList, i))) {
            IsFactor[i] = 1;
        } else {
            IsFactor[i] = 0;
        }

        sumLength += Rf_length(VECTOR_ELT(RList, i));
    }

    std::vector<int> primes;
    getAtLeastNPrimes(primes, sumLength);

    int numFactorVec = std::accumulate(IsFactor.cbegin(),
                                       IsFactor.cend(), 0);
    const int IsRep = CleanConvert::convertFlag(RIsRep, "IsRep");

    // All duplicates have been removed from RList via
    // lapply(RList, function(x) sort(unique(x)))
    std::vector<std::vector<int>> myVec(nCols);
    std::unordered_map<std::string, int> mapIndex;
    std::vector<int> typeCheck(5, 0);

    SEXP charVec     = PROTECT(Rf_allocVector(STRSXP, sumLength));
    SEXP dblSexpVec  = PROTECT(Rf_allocVector(REALSXP, sumLength));
    SEXP intSexpVec  = PROTECT(Rf_allocVector(INTSXP, sumLength));
    SEXP boolSexpVec = PROTECT(Rf_allocVector(LGLSXP, sumLength));

    double* dblVec = REAL(dblSexpVec);
    int* intVec    = INTEGER(intSexpVec);
    int* boolVec   = INTEGER(boolSexpVec);

    std::vector<std::vector<int>> facList(numFactorVec,
                                          std::vector<int>(sumLength, 0));

    for (int i = 0, total = 0, myIndex = 0, facInd = 0; i < nCols; ++i) {
        rcppType myType;
        std::vector<std::string> tempVec;
        convertToString(tempVec, VECTOR_ELT(RList, i), myType, IsFactor[i]);
        std::size_t j = 0;

        for (const auto &v: tempVec) {
            if (mapIndex.empty()) {
                mapIndex.insert({v, total});
            } else {
                auto mapCheck = mapIndex.find(v);

                if (mapCheck == mapIndex.end()) {
                    ++total;
                    mapIndex.insert({v, total});
                    myIndex = total;
                } else {
                    myIndex = mapCheck->second;
                }
            }

            myVec[i].push_back(myIndex);

            switch(myType) {
                case tInt: {
                    intVec[myIndex] = INTEGER(VECTOR_ELT(RList, i))[j];
                    typeCheck[tInt] = 1;
                    break;
                } case tFac: {
                    facList[facInd][myIndex] = INTEGER(VECTOR_ELT(RList, i))[j];
                    typeCheck[tFac] = 1;
                    break;
                } case tDbl: {
                    dblVec[myIndex] = REAL(VECTOR_ELT(RList, i))[j];
                    typeCheck[tDbl] = 1;
                    break;
                } case tStr: {
                    SET_STRING_ELT(charVec, myIndex, STRING_ELT(VECTOR_ELT(RList, i), j));
                    typeCheck[tStr] = 1;
                    break;
                } case tLog: {
                    boolVec[myIndex] = INTEGER(VECTOR_ELT(RList, i))[j];
                    typeCheck[tLog] = 1;
                    break;
                }
            }

            ++j;
        }

        facInd += IsFactor[i];
    }

    int mySum = std::accumulate(typeCheck.cbegin(), typeCheck.cend(), 0);

    // We need to check to see if there is overlap in factor levels
    if (typeCheck[tFac] && mySum == 1) {
        std::vector<std::string> testLevels;

        for (int i = 0; i < nCols; ++i) {
            if (IsFactor[i]) {
                SEXP facVec = PROTECT(Rf_getAttrib(VECTOR_ELT(RList, i),
                                                   R_LevelsSymbol));
                std::vector<std::string> strVec;

                int len_comp = Rf_length(facVec);

                for (int i = 0; i < len_comp; ++i) {
                    const std::string temp(CHAR(STRING_ELT(facVec, i)));
                    strVec.push_back(temp);
                }

                UNPROTECT(1);

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
    std::vector<int> lenGrps;
    std::vector<int> lastCol;
    std::vector<int> cartCombs;
    comboGrid(cartCombs, lastCol, lenGrps, myVec, primes, IsRep);

    const int nRows = lastCol.size();
    SEXP res =  PROTECT(GlueComboCart(cartCombs, lastCol, lenGrps, facList,
                                      typeCheck, IsFactor, RList, intVec,
                                      boolVec, dblVec, charVec, nRows, nCols, IsDF));
    UNPROTECT(5);
    return res;
}
