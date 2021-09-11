#include "ImportExportMPZ.h"
#include "CleanConvert.h"

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
        CleanConvert::convertVector(RindexVec, mySample,
                                    VecType::Numeric,
                                    "indexVec", false);
        sampSize = mySample.size();
        
        double myMax = *std::max_element(mySample.cbegin(), mySample.cend());
        
        if (myMax > computedRows) {
            Rf_error("One or more of the requested values exceeds"
                     " the maximum number of possible results");
        }

        if (sampSize > std::numeric_limits<int>::max()) {
            Rf_error("The number of rows cannot exceed 2^31 - 1");
        }

        // Get zero base index
        for (auto &s: mySample) {
            --s;
        }
    }
}

void SetIndexVecMpz(SEXP RindexVec, mpz_t *myVec,
                    std::size_t sampSize, mpz_t &computedRowsMpz) {
    
    createMPZArray(RindexVec, myVec, sampSize, "sampleVec");
    
    // get zero base
    for (std::size_t i = 0; i < sampSize; ++i)
        mpz_sub_ui(myVec[i], myVec[i], 1);
    
    mpz_t maxGmp;
    mpz_init(maxGmp);
    mpz_set(maxGmp, myVec[0]);
    
    for (std::size_t i = 1; i < sampSize; ++i)
        if (mpz_cmp(myVec[i], maxGmp) > 0)
            mpz_set(maxGmp, myVec[i]);
        
    if (mpz_cmp(maxGmp, computedRowsMpz) >= 0) {
        Rf_error("One or more of the requested values in sampleVec "
                 "exceeds the maximum number of possible results");
    }
}

void increment(bool IsGmp, mpz_t &mpzIndex, double &dblIndex) {
    if (IsGmp) {
        mpz_add_ui(mpzIndex, mpzIndex, 1u);
    } else {
        ++dblIndex;
    }
}

void increment(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int nRows) {
    if (IsGmp) {
        mpz_add_ui(mpzIndex, mpzIndex, nRows);
    } else {
        dblIndex += nRows;
    }
}

void decrement(bool IsGmp, mpz_t &mpzIndex, double &dblIndex) {
    if (IsGmp) {
        mpz_sub_ui(mpzIndex, mpzIndex, 1u);
    } else {
        --dblIndex;
    }
}

void decrement(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int nRows) {
    if (IsGmp) {
        mpz_sub_ui(mpzIndex, mpzIndex, nRows);
    } else {
        dblIndex -= nRows;
    }
}

bool CheckEqSi(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int si) {
    if (IsGmp) {
        return mpz_cmp_si(mpzIndex, si) == 0;
    } else {
        return dblIndex == si;
    }
}

bool CheckIndLT(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                mpz_t computedRowsMpz, double computedRows, bool eq = false) {
    if (eq) {
        if (IsGmp) {
            return mpz_cmp(mpzIndex, computedRowsMpz) <= 0;
        } else {
            return dblIndex <= computedRows;
        }
    } else {
        if (IsGmp) {
            return mpz_cmp(mpzIndex, computedRowsMpz) < 0;
        } else {
            return dblIndex < computedRows;
        }
    }
}

bool CheckEqInd(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                mpz_t computedRowsMpz, double computedRows) {
    if (IsGmp) {
        return mpz_cmp(mpzIndex, computedRowsMpz) == 0;
    } else {
        return dblIndex == computedRows;
    }
}

bool CheckIndGrT(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                 mpz_t computedRowsMpz, double computedRows) {
    if (IsGmp) {
        return mpz_cmp(mpzIndex, computedRowsMpz) > 0;
    } else {
        return dblIndex > computedRows;
    }
}

bool CheckGrTSi(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int si) {
    if (IsGmp) {
        return mpz_cmp_si(mpzIndex, si) > 0;
    } else {
        return dblIndex > si;
    }
}

template <typename T>
void UpdateExact(T* mat, T* yPt, T* xPt, std::vector<int> &z,
                 int lastRow, int nRows, int m, int n1) {
    
    for (int i = 0; i < m; ++i) {
        yPt[i] = mat[lastRow + i * nRows];
    }
    
    for (int j = 0; j < m; ++j) {
        int ind = 0;
        
        while (ind < n1 && xPt[ind] != yPt[j]) {
            ++ind;
        }
        
        z[j] = ind;
    }
}

void zUpdateIndex(std::vector<int> &z, SEXP v, SEXP mat, int m, int nRows) {
    
    constexpr double myTolerance = 8 * std::numeric_limits<double>::epsilon();
    const int n1 = Rf_length(v) - 1;
    const int lastRow = nRows - 1;
    
    SEXP vCopy = PROTECT(Rf_duplicate(v));
    z.assign(m, 0);
    
    switch (TYPEOF(v)) {
        case LGLSXP: {
            SEXP yBool = PROTECT(Rf_allocVector(LGLSXP, m));
            int* matBool = INTEGER(mat);
            int* xBoolPt = INTEGER(vCopy);
            int* yBoolPt = INTEGER(yBool);
            UpdateExact(matBool, yBoolPt, xBoolPt, z, lastRow, nRows, m, n1);
            break;
        } case INTSXP: {
            SEXP yInt = PROTECT(Rf_allocVector(INTSXP, m));
            int* matInt = INTEGER(mat);
            int* xIntPt = INTEGER(vCopy);
            int* yIntPt = INTEGER(yInt);
            UpdateExact(matInt, yIntPt, xIntPt, z, lastRow, nRows, m, n1);
            break;
        } case REALSXP: {
            SEXP yNum = PROTECT(Rf_allocVector(REALSXP, m));
            double* matNum = REAL(mat);
            double* xNumPt = REAL(vCopy);
            double* yNumPt = REAL(yNum);
            
            for (int i = 0; i < m; ++i) {
                yNumPt[i] = matNum[lastRow + i * nRows];
            }

            for (int j = 0; j < m; ++j) {
                int ind = 0;
                
                while (ind < n1 &&
                       std::abs(xNumPt[ind] - yNumPt[j]) > myTolerance) {

                    ++ind;
                }

                z[j] = ind;
            }
            
            break;
        } case STRSXP: {
            SEXP yChar = PROTECT(Rf_allocVector(STRSXP, m));
            
            for (int i = 0; i < m; ++i) {
                SET_STRING_ELT(yChar, i,
                               STRING_ELT(mat, lastRow + i * nRows));
            }

            for (int j = 0; j < m; ++j) {
                int ind = 0;
                
                while (ind < n1 &&
                       STRING_ELT(vCopy, ind) != STRING_ELT(yChar, j)) {
                    ++ind;
                }
                
                z[j] = ind;
            }

            break;
        } case CPLXSXP: {
            SEXP yCmplx = PROTECT(Rf_allocVector(CPLXSXP, m));
            Rcomplex* matCmplx = COMPLEX(mat);
            Rcomplex* xCmplxPt = COMPLEX(vCopy);
            Rcomplex* yCmplxPt = COMPLEX(yCmplx);

            for (int i = 0; i < m; ++i) {
                yCmplxPt[i] = matCmplx[lastRow + i * nRows];
            }

            for (int j = 0; j < m; ++j) {
                int ind = 0;
                bool bTestImg = std::abs(xCmplxPt[ind].i - yCmplxPt[j].i) > myTolerance;
                bool bTestReal = std::abs(xCmplxPt[ind].r - yCmplxPt[j].r) > myTolerance;
                
                while (ind < n1 && (bTestImg || bTestReal)) {
                    ++ind;
                    bTestImg = std::abs(xCmplxPt[ind].i - yCmplxPt[j].i) > myTolerance;
                    bTestReal = std::abs(xCmplxPt[ind].r - yCmplxPt[j].r) > myTolerance;
                }
                
                z[j] = ind;
            }
            
            break;
        } case RAWSXP: {
            SEXP yRaw = PROTECT(Rf_allocVector(RAWSXP, m));
            Rbyte* matRaw = RAW(mat);
            Rbyte* xRawPt = RAW(vCopy);
            Rbyte* yRawPt = RAW(yRaw);
            UpdateExact(matRaw, yRawPt, xRawPt, z, lastRow, nRows, m, n1);
            break;
        } default:{
            Rf_error("Only atomic types are supported for v");
        }
    }
    
    UNPROTECT(2);
}
