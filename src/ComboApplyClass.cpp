#include "ClassUtils/ComboApplyClass.h"

SEXP ComboApply::VecApplyReturn() {
    
    SEXP vectorPass = PROTECT(Rf_allocVector(RTYPE, m));
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
    
    switch (RTYPE) {
        case LGLSXP: {
            int* ptrOut = INTEGER(vectorPass);
            int* ptrIn  = INTEGER(sexpVec);
            
            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }
            
            break;
        } case INTSXP: {
            int* ptrOut = INTEGER(vectorPass);
            int* ptrIn  = INTEGER(sexpVec);
            
            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }
            
            if (IsFactor) {
                Rf_setAttrib(vectorPass, R_ClassSymbol, myClass);
                Rf_setAttrib(vectorPass, R_LevelsSymbol, myLevels);
            }
            
            break;
        } case STRSXP: {
            for (int j = 0; j < m; ++j) {
            SET_STRING_ELT(vectorPass, j, STRING_ELT(sexpVec, z[j]));
        }
            
            break;
        } case CPLXSXP: {
            Rcomplex* ptrOut = COMPLEX(vectorPass);
            Rcomplex* ptrIn  = COMPLEX(sexpVec);
            
            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }
            
            break;
        } case RAWSXP: {
            Rbyte* ptrOut = RAW(vectorPass);
            Rbyte* ptrIn  = RAW(sexpVec);
            
            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }
            
            break;
        } default: {
            double* ptrOut = REAL(vectorPass);
            double* ptrIn  = REAL(sexpVec);
            
            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }
            
            break;
        }
    }
    
    SETCADR(sexpFun, vectorPass);
    SEXP res = Rf_eval(sexpFun, rho);
    UNPROTECT(2);
    return res;
}

SEXP ComboApply::ApplyForward(int nRows) {
    return GetCombPermApply(sexpVec, vNum, vInt, n, m, IsComb,
                            IsRep, IsMult, freqs, z, myReps,
                            myType, nRows, stdFun, rho, RFunVal);
}

SEXP ComboApply::ApplyReverse(int nRows) {
    return GetPrevCombPermApply(sexpVec, vNum, vInt, myReps, freqs, z,
                                prevIter, n, m, IsComb, IsMult, nRows,
                                myType, stdFun, rho, RFunVal);
}

// The bVec Vector represents IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
ComboApply::ComboApply(SEXP Rv, int Rm, SEXP RcompRows,
                       const std::vector<int> &bVec,
                       const std::vector<int> &Rreps,
                       const std::vector<int> &Rfreqs,
                       const std::vector<int> &RvInt,
                       const std::vector<double> &RvNum,
                       VecType typePass, int RmaxThreads,
                       SEXP RnumThreads, SEXP RstdFun,
                       SEXP Rrho, SEXP R_RFunVal)
    : Combo(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum,
            typePass, RmaxThreads, RnumThreads), rho(Rrho),
            stdFun(RstdFun), RFunVal(R_RFunVal) {}

SEXP ComboApply::nextComb() {
    
    if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        return VecApplyReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex,
                          computedRowsMpz[0], computedRows)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        nextIter(freqs, z, n1, m1);
        return VecApplyReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex,
                          computedRowsMpz[0], computedRows)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::prevComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return VecApplyReturn();
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        prevIter(freqs, z, n1, m1);
        return VecApplyReturn();
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        const std::string message = "Iterator Initialized. To see the first"
                                    " result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::nextNumCombs(SEXP RNum) {
    
    int num;
    CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");
    
    if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex,
                   computedRowsMpz[0], computedRows)) {
        
        int nRows = 0;
        int numIncrement = 0;
        
        if (IsGmp) {
            mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numIncrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = computedRows - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
            nextIter(freqs, z, n1, m1);
        }
        
        increment(IsGmp, mpzIndex[0], dblIndex, numIncrement);
        SEXP res = PROTECT(ApplyForward(nRows));

        if (IsGmp) {
            mpz_sub_ui(mpzTemp, mpzIndex[0], 1u);
        } else {
            dblTemp = dblIndex - 1;
        }

        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        const std::string message = "No more results. To see the last result"
                                    ", use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::prevNumCombs(SEXP RNum) {
    
    int num;
    CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");
    
    if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        int nRows = 0;
        int numDecrement = 0;
        
        if (IsGmp) {
            mpz_sub_ui(mpzTemp, mpzIndex[0], 1u);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numDecrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = dblIndex - 1;
            nRows = num > dblTemp ? dblTemp : num;
            numDecrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex,
                       computedRowsMpz[0], computedRows, true)) {
            prevIter(freqs, z, n1, m1);
        }

        decrement(IsGmp, mpzIndex[0], dblIndex, numDecrement);
        return ApplyReverse(nRows);
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        const std::string message = "No more results. To see the last result"
                                    ", use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::nextGather() {
    
    if (IsGmp) {
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
        
        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = computedRows - dblIndex;
        
        if (dblTemp > std::numeric_limits<int>::max()) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }
    
    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;
    
    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0))
            nextIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], computedRowsMpz[0], 1u);
        } else {
            dblIndex = computedRows + 1;
        }
        
        SEXP res = PROTECT(ApplyForward(nRows));
        
        if (IsGmp) {
            mpz_sub_ui(mpzTemp, computedRowsMpz[0], 1u);
        } else {
            dblTemp = computedRows - 1;
        }
        
        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::prevGather() {
    
    if (IsGmp) {
        mpz_sub_ui(mpzTemp, mpzIndex[0], 1);
        
        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = dblIndex - 1;
        
        if (dblTemp > std::numeric_limits<int>::max()) {
            Rf_error("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }
    
    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;
    
    if (nRows) {
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_set_si(mpzIndex[0], 0u);
        } else {
            dblIndex = 0;
        }
        
        return ApplyReverse(nRows);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboApply::currComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex,
                    computedRowsMpz[0], computedRows)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        return VecApplyReturn();
    } else {
        const std::string message = "Iterator Initialized. To see the first "
                                    "result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    }
}

SEXP ComboApply::randomAccess(SEXP RindexVec) {
    
    std::size_t sampSize;
    std::vector<double> mySample;
    SetIndexVec(RindexVec, mySample, sampSize, IsGmp, computedRows);
    
    const std::size_t bigSampSize = (IsGmp) ? sampSize : 1;
    auto mpzVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i) {
        mpz_init(mpzVec[i]);
    }
    
    if (IsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec.get(), sampSize, computedRowsMpz[0]);
    }
    
    if (sampSize > 1) {
        return SampleCombPermApply(
            sexpVec, vInt, vNum, mySample, mpzVec.get(), myReps, stdFun,
            rho, RFunVal, nthResFun, myType, n, m, sampSize, false, IsGmp
        );
    } else {
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], mpzVec[0], 1u);
            mpz_set(mpzTemp, mpzVec[0]);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp = mySample.front();
        }
        
        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        return VecApplyReturn();
    }
}

SEXP ComboApply::front() {
    
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecApplyReturn();
}

SEXP ComboApply::back() {
    
    if (IsGmp) {
        mpz_set(mpzIndex[0], computedRowsMpz[0]);
        mpz_sub_ui(mpzTemp, computedRowsMpz[0], 1u);
    } else {
        dblIndex = computedRows;
        dblTemp = computedRows - 1;
    }
    
    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecApplyReturn();
}
