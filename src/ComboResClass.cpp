#include "ClassUtils/ComboResClass.h"

SEXP ComboRes::ApplyFun(SEXP mat) {

    if (Rf_isLogical(mat)) {
        return mat;
    }

    const int nRows = Rf_nrows(mat);
    SEXP res = PROTECT(Rf_allocMatrix(RTYPE, nRows, m + 1));

    if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);
        int* ptrIn  = INTEGER(mat);
        std::vector<int> vPass(m);

        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = ptrIn[i + nRows * j];
                ptrOut[i + nRows * j] = vPass[j];
            }

            ptrOut[i + nRows * m] = funInt(vPass, m);
        }
    } else {
        double* ptrOut = REAL(res);
        double* ptrIn  = REAL(mat);
        std::vector<double> vPass(m);

        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = ptrIn[i + nRows * j];
                ptrOut[i + nRows * j] = vPass[j];
            }

            ptrOut[i + nRows * m] = funDbl(vPass, m);
        }
    }

    UNPROTECT(1);
    return res;
}

SEXP ComboRes::VecReturn() {

    SEXP res = PROTECT(Rf_allocVector(RTYPE, m + 1));

    if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);
        std::vector<int> vPass(m);

        for (int j = 0; j < m; ++j) {
            vPass[j] = vInt[z[j]];
            ptrOut[j] = vPass[j];
        }

        ptrOut[m] = funInt(vPass, m);
    } else {
        double* ptrOut = REAL(res);
        std::vector<double> vPass(m);

        for (int j = 0; j < m; ++j) {
            vPass[j] = vNum[z[j]];
            ptrOut[j] = vPass[j];
        }

        ptrOut[m] = funDbl(vPass, m);
    }

    UNPROTECT(1);
    return res;
}

SEXP ComboRes::MatrixReturn(int nRows) {

    double userNum = 0;
    bool bSetNum = !numUnknown ||
        ctype == ConstraintType::SpecialCnstrnt;

    if (IsGmp) {
        mpz_add_ui(mpzTemp, mpzIndex, nRows);
    } else {
        dblTemp = dblIndex + nRows;
    }

    bLower = (IsGmp) ? mpz_cmp_si(mpzIndex, 0) > 0 : dblIndex > 0;
    bUpper = (IsGmp) ? mpz_cmp(mpzTemp, cnstrtCountMpz) < 0 :
      dblTemp < cnstrtCount;

    SetNumResults(IsGmp, bLower, bUpper, bSetNum, mpzTemp,
                  mpzIndex, dblIndex, dblTemp, cnstrtCount,
                  cnstrtCountMpz, nRows, userNum);

    int nThreads = 1;
    bool LocalPar = Parallel;
    const int limit = (part.isPart) ?
    ((part.ptype == PartitionType::RepCapped   ||
      part.ptype == PartitionType::DstctCapped ||
      part.ptype == PartitionType::DstctCappedMZ) ? 150000 : 40000) : 20000;

    SetThreads(LocalPar, maxThreads, nRows,
               myType, nThreads, sexpNThreads, limit);

    return GetConstraints(
        part, compVec, freqs, myReps, vNum, vInt, tarVals, tarIntVals, z,
        mainFun, funDbl, dblIndex, mpzIndex, userNum, ctype, myType, nThreads,
        nRows, n, strtLen, cap, m, IsComb, LocalPar, IsGmp, IsRep, IsMult,
        bUpper, KeepRes, numUnknown
    );
}

ComboRes::ComboRes(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
    std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
    std::vector<int> &RstartZ, const std::string &RmainFun,
    funcPtr<double> RfunDbl, ConstraintType Rctype, int RstrtLen,
    int Rcap, bool RKeepRes, bool RnumUnknown, double RcnstrtRows,
    mpz_t RcnstrtRowsMpz
) : Combo(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
          RmaxThreads, RnumThreads, Rparallel), cap(Rcap), strtLen(RstrtLen),
          KeepRes(RKeepRes), numUnknown(RnumUnknown), cnstrtCount(RcnstrtRows),
          tarIntVals(RtarIntVals), tarVals(RtarVals), ctype(Rctype),
          part(Rpart), mainFun(RmainFun), compVec(RcompVec),
          funDbl(RfunDbl), funInt(GetFuncPtr<int>(mainFun)) {

    z = RstartZ;
    mpz_init(cnstrtCountMpz);
    mpz_set(cnstrtCountMpz, RcnstrtRowsMpz);
    RTYPE = (myType == VecType::Integer) ? INTSXP : REALSXP;
}

void ComboRes::startOver() {
    Combo::startOver();
}

SEXP ComboRes::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return VecReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextIter(freqs, z, n1, m1);
        return VecReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex, dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboRes::prevComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    cnstrtCountMpz, cnstrtCount)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        return VecReturn();
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 1)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        prevIter(freqs, z, n1, m1);
        return VecReturn();
    } else if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        const std::string message = "Iterator Initialized. To see the first"
                                    " result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        decrement(IsGmp, mpzIndex, dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboRes::nextNumCombs(SEXP RNum) {

    int num;
    CleanConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {

        int nRows = 0;
        int numIncrement = 0;

        if (IsGmp) {
            mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numIncrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = cnstrtCount - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }

        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextIter(freqs, z, n1, m1);
        }

        SEXP res = PROTECT(MatrixReturn(nRows));
        increment(IsGmp, mpzIndex, dblIndex, numIncrement);
        zUpdateIndex(vNum, vInt, z, sexpVec, res, m, nRows);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last result"
                                    ", use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        increment(IsGmp, mpzIndex, dblIndex);
        return Rf_ScalarLogical(false);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboRes::prevNumCombs(SEXP RNum) {
    SEXP mat = PROTECT(Combo::prevNumCombs(RNum));
    SEXP res = PROTECT(ApplyFun(mat));
    UNPROTECT(2);
    return res;
}

SEXP ComboRes::nextGather() {

    if (IsGmp) {
        mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);

        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rf_error("The number of requested rows is greater than ",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = cnstrtCount - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            Rf_error("The number of requested rows is greater than ",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextIter(freqs, z, n1, m1);
        }

        SEXP res = PROTECT(MatrixReturn(nRows));

        if (IsGmp) {
            mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
        } else {
            dblIndex = cnstrtCount + 1;
        }

        zUpdateIndex(vNum, vInt, z, sexpVec, res, m, nRows);
        if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP ComboRes::prevGather() {
    SEXP mat = PROTECT(Combo::prevGather());
    SEXP res = PROTECT(ApplyFun(mat));
    UNPROTECT(2);
    return res;
}

SEXP ComboRes::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    cnstrtCountMpz, cnstrtCount)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return VecReturn();
    } else {
        const std::string message = "Iterator Initialized. To see the first "
                                    "result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    }
}

SEXP ComboRes::randomAccess(SEXP RindexVec) {
    SEXP samp = PROTECT(Combo::randomAccess(RindexVec));
    SEXP res  = PROTECT(Rf_isMatrix(samp) ? ApplyFun(samp) : VecReturn());
    UNPROTECT(2);
    return res;
}

SEXP ComboRes::front() {

    if (IsGmp) {
        mpz_set_ui(mpzIndex, 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecReturn();
}

SEXP ComboRes::back() {

    if (IsGmp) {
        mpz_set(mpzIndex, cnstrtCountMpz);
        mpz_sub_ui(mpzTemp, cnstrtCountMpz, 1u);
    } else {
        dblIndex = cnstrtCount;
        dblTemp = cnstrtCount - 1;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecReturn();
}
