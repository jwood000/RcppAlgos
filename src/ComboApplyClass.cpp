#include "ClassUtils/ComboApplyClass.h"

SEXP ComboApply::VecApplyReturn() {

    cpp11::sexp vectorPass = Rf_allocVector(RTYPE, m);
    cpp11::sexp sexpFun    = Rf_lang2(stdFun, R_NilValue);

    switch (RTYPE) {
        case LGLSXP:
        case INTSXP: {
            int* ptrOut = INTEGER(vectorPass);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = vInt[z[j]];
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

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = vNum[z[j]];
            }

            break;
        }
    }

    SETCADR(sexpFun, Rf_duplicate(vectorPass));
    cpp11::sexp res = Rf_eval(sexpFun, rho);
    return res;
}

SEXP ComboApply::ApplyForward(int nRows) {
    return GetCombPermApply(sexpVec, vNum, vInt, n, m, IsComb,
                            IsRep, IsMult, freqs, z, myReps,
                            myType, nRows, stdFun, rho, RFunVal);
}

SEXP ComboApply::ApplyReverse(int nRows) {
    return GetPrevCombPermApply(sexpVec, vNum, vInt, myReps, freqs, z,
                                prevComb, n, m, IsComb, IsMult, nRows,
                                myType, stdFun, rho, RFunVal);
}

ComboApply::ComboApply(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal
) : Combo(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
          RmaxThreads, RnumThreads, Rparallel), rho(Rrho), stdFun(RstdFun),
          RFunVal(R_RFunVal) {}

void ComboApply::startOver() {
    Combo::startOver();
}

SEXP ComboApply::nextIter() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0) &&
        CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return VecApplyReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextComb(freqs, z, n1, m1);
        return VecApplyReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::prevIter() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        return VecApplyReturn();
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 1)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        prevComb(freqs, z, n1, m1);
        return VecApplyReturn();
    } else if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::nextNumIters(SEXP RNum) {

    int num;
    CppConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {

        int nRows = 0;
        int numIncrement = 0;

        if (IsGmp) {
            mpzTemp = computedRowsMpz - mpzIndex;
            nRows = cmp(mpzTemp, num) < 0 ? mpzTemp.get_si() : num;
            numIncrement = cmp(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = computedRows - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }

        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextComb(freqs, z, n1, m1);
        }

        increment(IsGmp, mpzIndex, dblIndex, numIncrement);
        cpp11::sexp res = ApplyForward(nRows);

        if (IsGmp) {
            mpzTemp = mpzIndex - 1;
        } else {
            dblTemp = dblIndex - 1;
        }

        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        return res;
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::prevNumIters(SEXP RNum) {

    int num;
    CppConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 1)) {
        int nRows = 0;
        int numDecrement = 0;

        if (IsGmp) {
            mpzTemp = mpzIndex - 1;
            nRows = cmp(mpzTemp, num) < 0 ? mpzTemp.get_si() : num;
            numDecrement = cmp(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = dblIndex - 1;
            nRows = num > dblTemp ? dblTemp : num;
            numDecrement = num > dblTemp ? (nRows + 1) : nRows;
        }

        if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                       computedRowsMpz, computedRows, true)) {
            prevComb(freqs, z, n1, m1);
        }

        decrement(IsGmp, mpzIndex, dblIndex, numDecrement);
        return ApplyReverse(nRows);
    } else if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::nextGather() {

    if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        return ToSeeLast();
    }

    if (IsGmp) {
        mpzTemp = computedRowsMpz - mpzIndex;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = computedRows - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)){
            nextComb(freqs, z, n1, m1);
        }

        if (IsGmp) {
            mpzIndex = computedRowsMpz + 1;
        } else {
            dblIndex = computedRows + 1;
        }

        cpp11::sexp res = ApplyForward(nRows);

        if (IsGmp) {
            mpzTemp = computedRowsMpz - 1;
        } else {
            dblTemp = computedRows - 1;
        }

        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        return res;
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::prevGather() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    }

    if (IsGmp) {
        mpzTemp = mpzIndex - 1;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = dblIndex - 1;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of requested rows is greater than ",
                     std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                       computedRowsMpz, computedRows, true)) {
            prevComb(freqs, z, n1, m1);
        }

        if (IsGmp) {
            mpzIndex = 0;
        } else {
            dblIndex = 0;
        }

        return ApplyReverse(nRows);
    } else {
        return R_NilValue;
    }
}

SEXP ComboApply::currIter() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return VecApplyReturn();
    } else {
        return ToSeeFirst(false);
    }
}

SEXP ComboApply::randomAccess(SEXP RindexVec) {

    std::size_t sampSize;
    std::vector<double> mySample;
    SetIndexVec(RindexVec, mySample, sampSize, IsGmp, computedRows);

    const std::size_t bigSampSize = IsGmp ? sampSize : 1;
    std::vector<mpz_class> mpzVec(bigSampSize);

    if (IsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec, sampSize, computedRowsMpz);
    }

    if (sampSize > 1) {
        return SampleCombPermApply(
            sexpVec, vInt, vNum, mySample, mpzVec, myReps, stdFun,
            rho, RFunVal, nthResFun, myType, n, m, sampSize, false, IsGmp
        );
    } else {
        if (IsGmp) {
            mpzIndex = mpzVec.front() + 1;
            mpzTemp  = mpzVec.front();
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp  = mySample.front();
        }

        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        return VecApplyReturn();
    }
}

SEXP ComboApply::front() {

    if (IsGmp) {
        mpzIndex = 1;
        mpzTemp  = 0;
    } else {
        dblIndex = 1;
        dblTemp  = 0;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecApplyReturn();
}

SEXP ComboApply::back() {

    if (IsGmp) {
        mpzIndex = computedRowsMpz;
        mpzTemp  = computedRowsMpz - 1;
    } else {
        dblIndex = computedRows;
        dblTemp  = computedRows - 1;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecApplyReturn();
}
