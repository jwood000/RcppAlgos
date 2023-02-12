#include "Constraints/ConstraintsTypes.h"
#include "SetUpUtils.h"
#include <algorithm>   // std::max_element
#include <cmath>       // std::abs

void SetIndexVec(SEXP RindexVec, std::vector<double> &mySample,
                 std::size_t &sampSize, bool IsGmp, double computedRows) {

    if (IsGmp) {
        switch (TYPEOF(RindexVec)) {
            case RAWSXP: {
                const char* raw = (char*) RAW(RindexVec);
                sampSize = ((int*) raw)[0];
                break;
            } default: {
                sampSize = LENGTH(RindexVec);
            }
        }
    } else {
        CppConvert::convertVector(RindexVec, mySample, VecType::Numeric,
                                  "indexVec", false);
        sampSize = mySample.size();

        double myMax = *std::max_element(mySample.cbegin(), mySample.cend());

        if (myMax > computedRows) {
            cpp11::stop("One or more of the requested values exceeds"
                     " the maximum number of possible results");
        }

        if (sampSize > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of rows cannot exceed 2^31 - 1");
        }

        // Get zero base index
        for (auto &s: mySample) {
            --s;
        }
    }
}

void SetIndexVecMpz(SEXP RindexVec, std::vector<mpz_class> &myVec,
                    std::size_t sampSize, mpz_class computedRowsMpz) {

    CppConvert::convertMPZVector(RindexVec, myVec, sampSize, "sampleVec");

    // get zero base
    for (std::size_t i = 0; i < sampSize; ++i) {
        myVec[i]--;
    }

    mpz_class maxGmp(myVec[0]);

    for (std::size_t i = 1; i < sampSize; ++i) {
        if (cmp(myVec[i], maxGmp) > 0) {
            maxGmp = myVec[i];
        }
    }

    if (cmp(maxGmp, computedRowsMpz) >= 0) {
        cpp11::stop("One or more of the requested values in sampleVec "
                 "exceeds the maximum number of possible results");
    }
}

void increment(bool IsGmp, mpz_class &mpzIndex, double &dblIndex) {
    if (IsGmp) {
        ++mpzIndex;
    } else {
        ++dblIndex;
    }
}

void increment(bool IsGmp, mpz_class &mpzIndex, double &dblIndex, int nRows) {
    if (IsGmp) {
        mpzIndex += nRows;
    } else {
        dblIndex += nRows;
    }
}

void decrement(bool IsGmp, mpz_class &mpzIndex, double &dblIndex) {
    if (IsGmp) {
        --mpzIndex;
    } else {
        --dblIndex;
    }
}

void decrement(bool IsGmp, mpz_class &mpzIndex, double &dblIndex, int nRows) {
    if (IsGmp) {
        mpzIndex -= nRows;
    } else {
        dblIndex -= nRows;
    }
}

bool CheckEqSi(bool IsGmp, const mpz_class &mpzIndex,
               double dblIndex, int si) {
    if (IsGmp) {
        return cmp(mpzIndex, si) == 0;
    } else {
        return dblIndex == si;
    }
}

bool CheckIndLT(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                const mpz_class &computedRowsMpz, double computedRows,
                bool eq = false) {
    if (eq) {
        if (IsGmp) {
            return cmp(mpzIndex, computedRowsMpz) <= 0;
        } else {
            return dblIndex <= computedRows;
        }
    } else {
        if (IsGmp) {
            return cmp(mpzIndex, computedRowsMpz) < 0;
        } else {
            return dblIndex < computedRows;
        }
    }
}

bool CheckEqInd(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                const mpz_class &computedRowsMpz, double computedRows) {
    if (IsGmp) {
        return cmp(mpzIndex, computedRowsMpz) == 0;
    } else {
        return dblIndex == computedRows;
    }
}

bool CheckIndGrT(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                 const mpz_class &computedRowsMpz, double computedRows) {
    if (IsGmp) {
        return cmp(mpzIndex, computedRowsMpz) > 0;
    } else {
        return dblIndex > computedRows;
    }
}

bool CheckGrTSi(bool IsGmp, const mpz_class &mpzIndex,
                double dblIndex, int si) {
    if (IsGmp) {
        return cmp(mpzIndex, si) > 0;
    } else {
        return dblIndex > si;
    }
}

template <typename T>
void UpdateExact(T* mat, T* yPt, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t lastRow,
                 std::size_t nRows, std::size_t m, int n1, int numAdd = 0) {

    for (std::size_t j = 0; j < m; ++j) {
        yPt[j] = mat[lastRow + j * nRows];
    }

    for (std::size_t j = 0; j < m; ++j) {
        int ind = 0;

        while (ind < n1 && v[ind] != yPt[j]) {
            ++ind;
        }

        z[j] = ind + numAdd;
    }
}

void zUpdateIndex(const std::vector<double> &vNum,
                  const std::vector<int> &vInt, std::vector<int> &z,
                  SEXP v, SEXP mat, std::size_t m,
                  std::size_t nRows, bool bAddOne) {

    constexpr double myTolerance = 8 * std::numeric_limits<double>::epsilon();
    const int n1 = Rf_length(v) - 1;
    const std::size_t lastRow = nRows - 1;
    z.assign(static_cast<int>(m), 0);

    switch (TYPEOF(mat)) {
        case LGLSXP: {
            cpp11::sexp yBool = Rf_allocVector(LGLSXP, m);
            int* matBool = INTEGER(mat);
            int* yBoolPt = INTEGER(yBool);
            UpdateExact(matBool, yBoolPt, vInt, z, lastRow, nRows, m, n1);
            break;
        } case INTSXP: {
            const int numAdd = static_cast<int>(bAddOne);
            cpp11::sexp yInt = Rf_allocVector(INTSXP, m);
            int* matInt = INTEGER(mat);
            int* yIntPt = INTEGER(yInt);
            UpdateExact(matInt, yIntPt, vInt, z,
                        lastRow, nRows, m, n1, numAdd);
            break;
        } case REALSXP: {
            cpp11::sexp yNum = Rf_allocVector(REALSXP, m);
            double* matNum = REAL(mat);
            double* yNumPt = REAL(yNum);

            for (std::size_t i = 0; i < m; ++i) {
                yNumPt[i] = matNum[lastRow + i * nRows];
            }

            for (std::size_t j = 0; j < m; ++j) {
                int ind = 0;

                while (ind < n1 &&
                       std::abs(vNum[ind] - yNumPt[j]) > myTolerance) {

                    ++ind;
                }

                z[j] = ind;
            }

            break;
        } case STRSXP: {
            cpp11::sexp yChar = Rf_allocVector(STRSXP, m);

            for (std::size_t i = 0; i < m; ++i) {
                SET_STRING_ELT(yChar, i,
                               STRING_ELT(mat, lastRow + i * nRows));
            }

            for (std::size_t j = 0; j < m; ++j) {
                int ind = 0;

                while (ind < n1 &&
                       STRING_ELT(v, ind) != STRING_ELT(yChar, j)) {
                    ++ind;
                }

                z[j] = ind;
            }

            break;
        } case CPLXSXP: {
            cpp11::sexp yCmplx = Rf_allocVector(CPLXSXP, m);
            Rcomplex* matCmplx = COMPLEX(mat);
            Rcomplex* xCmplxPt = COMPLEX(v);
            Rcomplex* yCmplxPt = COMPLEX(yCmplx);

            for (std::size_t i = 0; i < m; ++i) {
                yCmplxPt[i] = matCmplx[lastRow + i * nRows];
            }

            for (std::size_t j = 0; j < m; ++j) {
                int ind = 0;
                bool bTestImg  = std::abs(xCmplxPt[ind].i - yCmplxPt[j].i) > myTolerance;
                bool bTestReal = std::abs(xCmplxPt[ind].r - yCmplxPt[j].r) > myTolerance;

                while (ind < n1 && (bTestImg || bTestReal)) {
                    ++ind;
                    bTestImg  = std::abs(xCmplxPt[ind].i - yCmplxPt[j].i) > myTolerance;
                    bTestReal = std::abs(xCmplxPt[ind].r - yCmplxPt[j].r) > myTolerance;
                }

                z[j] = ind;
            }

            break;
        } case RAWSXP: {
            cpp11::sexp yRaw = Rf_allocVector(RAWSXP, m);
            Rbyte* matRaw = RAW(mat);
            Rbyte* xRawPt = RAW(v);
            std::vector<Rbyte> stlRawVec(n1 + 1);

            for (int i = 0; i <= n1; ++i) {
                stlRawVec[i] = xRawPt[i];
            }

            Rbyte* yRawPt = RAW(yRaw);
            UpdateExact(matRaw, yRawPt, stlRawVec, z, lastRow, nRows, m, n1);
            break;
        } default:{
            cpp11::stop("Only atomic types are supported for v");
        }
    }
}
