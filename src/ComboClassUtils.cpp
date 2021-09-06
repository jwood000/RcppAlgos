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
