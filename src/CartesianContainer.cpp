#include "cpp11/logicals.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"

#include "NumbersUtils/Eratosthenes.h"
#include "ComboCartesian.h"
#include "CppConvert.h"
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
                cpp11::integers facVec(ListElement);
                cpp11::strings myLevels(facVec.attr("levels"));
                std::vector<std::string> strVec;

                for (int i = 0, faclen = Rf_length(myLevels); i < faclen; ++i) {
                    strVec.push_back(myLevels[i]);
                }

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
                cpp11::integers intVec(ListElement);
                typePass = tInt;

                for (auto v: intVec) {
                    tempVec.push_back(std::to_string(v));
                }
            }

            break;
        } case LGLSXP: {
            cpp11::logicals boolVec(ListElement);
            typePass = tLog;

            for (auto v: boolVec) {
                std::string r = std::to_string(v ? 1 : 0);
                tempVec.push_back(r + "log");
            }

            break;
        } case REALSXP: {
            cpp11::doubles dblVec(ListElement);
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

void GetCharOutput(cpp11::writable::strings_matrix<> &res,
                   const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const cpp11::strings &charVec,
                   int nCols, int nRows) {

    for (int i = 0, n1 = nCols - 1, row = 0, m_idx = 0,
         baseSize = lenGrps.size(); i < baseSize; ++i, m_idx = row) {

        for (int j = 0, c_idx = i * n1, grpSize = lenGrps[i];
             j < n1; ++j, ++c_idx, m_idx += nRows) {

            // Benchmarks show that this set up is about 30%
            // faster than using cpp11::sexp
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
}

template <typename matType, typename T>
void GetPureOutput(matType &mat,
                   const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const std::vector<T> &standardVec,
                   int nCols, int nRows) {

    for (int i = 0, n1 = nCols - 1, row = 0, m_idx = 0,
         baseSize = lenGrps.size(); i < baseSize; ++i) {

        for (int j = 0, c_idx = i * n1, grpSize = lenGrps[i];
             j < n1; ++j, ++c_idx) {

            auto&& comb = standardVec[cartCombs[c_idx]];

            for (int k = 0; k < grpSize; ++k) {
                mat(m_idx + k, j) = comb;
            }
        }

        for (int k = 0; k < lenGrps[i]; ++k, ++row) {
            mat(m_idx + k, n1) = standardVec[lastCol[row]];
        }

        m_idx += lenGrps[i];
    }
}

template <typename cpp11Type, typename T>
void PoulateColumn(const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const std::vector<T> &poolVec,
                   cpp11Type &res, int nCols, int nRows, int i) {

    if (i < (nCols - 1)) {
        for (int j = 0, n1 = nCols - 1, row = i, idx = 0,
             baseSize = lenGrps.size(); idx < baseSize;
             ++idx, row += n1) {

            auto&& comb = poolVec[cartCombs[row]];

            for (int k = 0; k < lenGrps[idx]; ++k, ++j) {
                res[j] = comb;
            }
        }
    } else {
        for (int j = 0; j < nRows; ++j) {
            res[j] = poolVec[lastCol[j]];
        }
    }
}

SEXP GlueComboCart(const std::vector<int> &cartCombs,
                   const std::vector<int> &lastCol,
                   const std::vector<int> &lenGrps,
                   const std::vector<std::vector<int>> &facList,
                   const std::vector<int> &typeCheck,
                   const std::vector<int> &IsFactor,
                   const cpp11::list &RList,
                   const std::vector<int> &intVec,
                   const std::vector<double> &dblVec,
                   const std::vector<int> &boolVec,
                   const cpp11::strings &charVec,
                   int nRows, int nCols, bool IsDF) {

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);

        for (int i = 0, facInd = 0; i < nCols; ++i) {
            switch (TYPEOF(RList[i])) {
                case INTSXP: {
                    cpp11::writable::integers intSexpVec(nRows);

                    if (IsFactor[i]) {
                        const std::vector<int> facVec(
                            facList[facInd].begin(), facList[facInd].end()
                        );

                        PoulateColumn(cartCombs, lastCol, lenGrps,
                                      facVec, intSexpVec, nCols, nRows, i);
                        SetFactorClass(intSexpVec, RList[i]);
                        ++facInd;
                    } else {
                        PoulateColumn(cartCombs, lastCol, lenGrps,
                                      intVec, intSexpVec, nCols, nRows, i);
                    }

                    DataFrame[i] = intSexpVec;
                    break;
                } case LGLSXP: {
                    cpp11::writable::logicals boolSexpVec(nRows);
                    PoulateColumn(cartCombs, lastCol, lenGrps,
                                  boolVec, boolSexpVec, nCols, nRows, i);
                    DataFrame[i] = boolSexpVec;
                    break;
                } case REALSXP: {
                    cpp11::writable::doubles dblSexpVec(nRows);
                    PoulateColumn(cartCombs, lastCol, lenGrps,
                                  dblVec, dblSexpVec, nCols, nRows, i);
                    DataFrame[i] = dblSexpVec;
                    break;
                } case STRSXP: {
                    cpp11::writable::strings sexpVec(nRows);

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

                    DataFrame[i] = sexpVec;
                    break;
                }
            }
        }

        DataFrame.attr("row.names") = {NA_INTEGER, -nRows};
        DataFrame.names() = RList.names();
        DataFrame.attr("class") = "data.frame";
        return DataFrame;
    } else {
        if (typeCheck[tInt]) {
            cpp11::writable::integers_matrix<> intMat(nRows, nCols);
            GetPureOutput(intMat, cartCombs, lastCol,
                          lenGrps, intVec, nCols, nRows);
            return intMat;
        } else if (typeCheck[tFac]) {
            cpp11::writable::integers_matrix<> intMat(nRows, nCols);
            GetPureOutput(intMat, cartCombs, lastCol,
                          lenGrps, intVec, nCols, nRows);
            SetFactorClass(intMat, RList[0]);
            return intMat;
        } else if (typeCheck[tLog]) {
            cpp11::writable::logicals_matrix<> boolMat(nRows, nCols);
            GetPureOutput(boolMat, cartCombs, lastCol,
                          lenGrps, boolVec, nCols, nRows);
            return boolMat;
        } else if (typeCheck[tDbl]) {
            cpp11::writable::doubles_matrix<> dblMat(nRows, nCols);
            GetPureOutput(dblMat, cartCombs, lastCol,
                          lenGrps, dblVec, nCols, nRows);
            return dblMat;
        } else {
            cpp11::writable::strings_matrix<> res(nRows, nCols);
            GetCharOutput(res, cartCombs, lastCol,
                          lenGrps, charVec, nCols, nRows);
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

[[cpp11::register]]
SEXP ComboGridCpp(cpp11::list RList, bool IsRep) {

    int sumLength = 0;
    const int nCols = Rf_length(RList);
    std::vector<int> IsFactor(nCols);

    for (int i = 0; i < nCols; ++i) {
        if (Rf_isFactor(RList[i])) {
            IsFactor[i] = 1;
        } else {
            IsFactor[i] = 0;
        }

        sumLength += Rf_length(RList[i]);
    }

    std::vector<int> primes;
    getAtLeastNPrimes(primes, sumLength);

    int numFactorVec = std::accumulate(IsFactor.cbegin(),
                                       IsFactor.cend(), 0);

    // All duplicates have been removed from RList via
    // lapply(RList, function(x) sort(unique(x)))
    std::vector<std::vector<int>> myVec(nCols);
    std::unordered_map<std::string, int> mapIndex;
    std::vector<int> typeCheck(5, 0);

    cpp11::writable::strings charVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    std::vector<std::vector<int>> facList(numFactorVec,
                                          std::vector<int>(sumLength, 0));

    for (int i = 0, total = 0, myIndex = 0, facInd = 0; i < nCols; ++i) {
        rcppType myType;
        std::vector<std::string> tempVec;
        convertToString(tempVec, RList[i], myType, IsFactor[i]);
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
                    intVec[myIndex] = INTEGER(RList[i])[j];
                    typeCheck[tInt] = 1;
                    break;
                } case tFac: {
                    facList[facInd][myIndex] = INTEGER(RList[i])[j];
                    typeCheck[tFac] = 1;
                    break;
                } case tDbl: {
                    dblVec[myIndex] = REAL(RList[i])[j];
                    typeCheck[tDbl] = 1;
                    break;
                } case tStr: {
                    SET_STRING_ELT(charVec, myIndex, STRING_ELT(RList[i], j));
                    typeCheck[tStr] = 1;
                    break;
                } case tLog: {
                    boolVec[myIndex] = LOGICAL(RList[i])[j];
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
    std::vector<int> lenGrps;
    std::vector<int> lastCol;
    std::vector<int> cartCombs;
    comboGrid(cartCombs, lastCol, lenGrps, myVec, primes, IsRep);

    const int nRows = lastCol.size();
    cpp11::sexp res = GlueComboCart(cartCombs, lastCol, lenGrps, facList,
                                    typeCheck, IsFactor, RList, intVec,
                                    dblVec, boolVec, charVec, nRows, nCols, IsDF);
    return res;
}
