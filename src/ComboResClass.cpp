#include "ClassUtils/ComboResClass.h"

SEXP ComboRes::ApplyFun(SEXP mat) {

    if (Rf_isLogical(mat)) {
        return mat;
    }

    const int nRows = Rf_nrows(mat);
    SEXP res = PROTECT(Rf_allocMatrix(RTYPE, nRows, nCols));

    if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);
        int* ptrIn  = INTEGER(mat);
        std::vector<int> vPass(width);

        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < width; ++j) {
                vPass[j] = ptrIn[i + nRows * j];
                ptrOut[i + nRows * j] = vPass[j];
            }

            ptrOut[i + nRows * width] = funInt(vPass, width);
        }
    } else {
        double* ptrOut = REAL(res);
        double* ptrIn  = REAL(mat);
        std::vector<double> vPass(width);

        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < width; ++j) {
                vPass[j] = ptrIn[i + nRows * j];
                ptrOut[i + nRows * j] = vPass[j];
            }

            ptrOut[i + nRows * width] = funDbl(vPass, width);
        }
    }

    UNPROTECT(1);
    return res;
}

SEXP ComboRes::VecReturn() {

    SEXP res = PROTECT(Rf_allocVector(RTYPE, nCols));

    if (ctype == ConstraintType::PartStandard) {
        int* ptrOut = INTEGER(res);

        for (int j = 0; j < width; ++j) {
          ptrOut[j] = z[j];
        }

        if (KeepRes) ptrOut[width] = part.target;
    } else if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);
        std::vector<int> vPass(width);

        for (int j = 0; j < width; ++j) {
            vPass[j] = vInt[z[j]];
            ptrOut[j] = vPass[j];
        }

        if (KeepRes) {
            if (part.isPart) {
                ptrOut[width] = part.target;
            } else {
                ptrOut[width] = funInt(vPass, width);
            }
        }
    } else {
        double* ptrOut = REAL(res);
        std::vector<double> vPass(width);

        for (int j = 0; j < width; ++j) {
            vPass[j] = vNum[z[j]];
            ptrOut[j] = vPass[j];
        }

        if (KeepRes) {
            if (part.isPart) {
                ptrOut[width] = part.target;
            } else {
                ptrOut[width] = funDbl(vPass, width);
            }
        }
    }

    UNPROTECT(1);
    return res;
}

SEXP ComboRes::MatrixReturn(int nRows) {

    dblTemp = 0;
    double userNum = nRows;
    mpz_set_ui(mpzTemp, 0u);

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
        mainFun, funTest, funDbl, dblTemp, mpzTemp, userNum, ctype, myType,
        nThreads, nRows, n, strtLen, cap, width, IsComb, LocalPar, IsGmp,
        IsRep, IsMult, bUpper, KeepRes, numUnknown
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
    const std::string &RFunTest, funcPtr<double> RfunDbl,
    ConstraintType Rctype, int RstrtLen, int Rcap, bool RKeepRes,
    bool RnumUnknown, double RcnstrtRows, mpz_t RcnstrtRowsMpz
) : Combo(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
          RmaxThreads, RnumThreads, Rparallel), cap(Rcap),
          width(Rpart.isPart ? Rpart.width : m), nCols(RKeepRes ? width + 1 :
                width), strtLen(RstrtLen),
          KeepRes(RKeepRes), numUnknown(RnumUnknown), cnstrtCount(RcnstrtRows),
          tarIntVals(RtarIntVals), tarVals(RtarVals), ctype(Rctype),
          part(Rpart), mainFun(RmainFun), funTest(RFunTest), compVec(RcompVec),
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
        return ToSeeLast();
    } else {
        return R_NilValue;
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
        return ToSeeFirst();
    } else {
        return R_NilValue;
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
            if (!nextIter(freqs, z, n1, m1)) {
                if (IsGmp) {
                    mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
                } else {
                    dblIndex = cnstrtCount + 1;
                }

                const std::string message = "No more results\n\n";
                Rprintf(message.c_str());
                return R_NilValue;
            }
        }

        SEXP res = PROTECT(MatrixReturn(nRows));
        increment(IsGmp, mpzIndex, dblIndex, numIncrement);
        zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows);
        if (!IsComb) TopOffPerm(z, myReps, n, width, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP ComboRes::prevNumCombs(SEXP RNum) {
    SEXP mat = PROTECT(Combo::prevNumCombs(RNum));
    SEXP res = PROTECT(ApplyFun(mat));
    UNPROTECT(2);
    return res;
}

SEXP ComboRes::nextGather() {

    if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast();
    }

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
            if (!nextIter(freqs, z, n1, m1)) {
                if (IsGmp) {
                    mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
                } else {
                    dblIndex = cnstrtCount + 1;
                }

                const std::string message = "No more results\n\n";
                Rprintf(message.c_str());
                return R_NilValue;
            }
        }

        SEXP res = PROTECT(MatrixReturn(nRows));

        if (IsGmp) {
            mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
        } else {
            dblIndex = cnstrtCount + 1;
        }

        zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows);
        if (!IsComb) TopOffPerm(z, myReps, n, width, IsRep, IsMult);
        UNPROTECT(1);
        return res;
    } else {
        return R_NilValue;
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
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return VecReturn();
    } else {
        return ToSeeFirst(false);
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

    z = nthResFun(n, width, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, width, IsRep, IsMult);
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

    z = nthResFun(n, width, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, width, IsRep, IsMult);
    return VecReturn();
}

SEXP ComboRes::summary() {
    SEXP res = PROTECT(Combo::summary());
    std::string desc(R_CHAR(STRING_ELT(VECTOR_ELT(res, 0), 0)));
    desc += " with " + mainFun + " applied to each result";
    SET_VECTOR_ELT(res, 0, Rf_mkString(desc.c_str()));
    UNPROTECT(1);
    return res;
}
