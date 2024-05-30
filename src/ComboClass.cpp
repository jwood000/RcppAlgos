#include "ClassUtils/ComboClass.h"

SEXP Combo::ToSeeLast(bool AdjustIdx) {
    std::string message = "No more results.";
    if (prevIterAvailable) {
        message += " To see the last result, use the prevIter method(s)\n\n";
    } else {
        message += "\n\n";
    }

    Rprintf("%s", message.c_str());
    if (AdjustIdx) increment(IsGmp, mpzIndex, dblIndex);
    return R_NilValue;
}

SEXP Combo::ToSeeFirst(bool AdjustIdx) {
    const std::string message = "Iterator Initialized. To see the first"
                                " result, use the nextIter method(s)\n\n";
    Rprintf("%s", message.c_str());
    if (AdjustIdx) decrement(IsGmp, mpzIndex, dblIndex);
    return R_NilValue;
}

SEXP Combo::BasicVecReturn() {

    cpp11::sexp res = Rf_allocVector(RTYPE, m);

    switch (RTYPE) {
        case LGLSXP:
        case INTSXP: {
            int* ptrOut = INTEGER(res);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = vInt[z[j]];
            }

            if (IsFactor) {
                Rf_setAttrib(res, R_ClassSymbol, myClass);
                Rf_setAttrib(res, R_LevelsSymbol, myLevels);
            }

            break;
        } case STRSXP: {
            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(res, j, STRING_ELT(sexpVec, z[j]));
            }

            break;
        } case CPLXSXP: {
            Rcomplex* ptrOut = COMPLEX(res);
            Rcomplex* ptrIn  = COMPLEX(sexpVec);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }

            break;
        } case RAWSXP: {
            Rbyte* ptrOut = RAW(res);
            Rbyte* ptrIn  = RAW(sexpVec);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }

            break;
        } default: {
            double* ptrOut = REAL(res);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = vNum[z[j]];
            }

            break;
        }
    }

    return res;
}

SEXP Combo::MatForward(int nRows, int numIncrement) {

    int nThreads = 1;
    bool LocalPar = Parallel;
    const int limit = 20000;

    SetThreads(LocalPar, maxThreads, nRows,
               myType, nThreads, sexpNThreads, limit);

    cpp11::sexp res = GetCombPerms(
        sexpVec, vNum, vInt, n, m, 0, true, IsComb, LocalPar, IsRep, IsMult,
        IsGmp, freqs, z, myReps, dblIndex, mpzIndex, nRows, nThreads, myType
    );

    zUpdateIndex(vNum, vInt, z, sexpVec, res, m, nRows);
    increment(IsGmp, mpzIndex, dblIndex, numIncrement);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return res;
}

SEXP Combo::MatReverse(int nRows) {
    return GetPrevCombPerms(sexpVec, vNum, vInt, myReps, freqs, z,
                            prevIter, n, m, IsComb, IsMult, nRows, myType);
}

// The bVec Vector represents IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
Combo::Combo(
    SEXP Rv, int Rm, SEXP RcompRow, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnThreads, bool Rparallel
) : n(Rf_length(Rv)), m(Rm), m1(Rm - 1), RTYPE(TYPEOF(Rv)),
    maxThreads(RmaxThreads), sexpVec(Rv), sexpNThreads(RnThreads),
    IsGmp(bVec[4]), IsFactor(bVec[0]), IsComb(bVec[1] && !bVec[6]),
    IsMult(bVec[2]), IsRep(bVec[3]), Parallel(Rparallel),
    computedRows(IsGmp ? 0 : Rf_asReal(RcompRow)), myType(typePass),
    vInt(RvInt), vNum(RvNum), freqs(Rfreqs), myReps(Rreps),
    n1(IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1)),
    myClass(bVec[0] ? Rf_getAttrib(Rv, R_ClassSymbol) :
                Rf_allocVector(STRSXP, 0)),
    myLevels(bVec[0] ? Rf_getAttrib(Rv, R_LevelsSymbol) : R_NilValue),
    nthResFun(GetNthResultFunc(bVec[1], bVec[2], bVec[3], bVec[4])),
    nextIter(GetNextIterPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
    prevIter(GetPrevIterPtr(bVec[1], bVec[2], bVec[3], bVec[5])) {

    z.resize(Rm);

    if (IsGmp) {
        CppConvert::convertMpzClass(RcompRow, computedRowsMpz,
                                    "computedRowsMpz");
    }

    dblIndex = 0;
    mpzIndex = 0;
    SetStartZ(myReps, freqs, z, IsComb, n, m, dblIndex,
              mpzIndex, IsRep, IsMult, IsGmp);
    prevIterAvailable = true;
}

void Combo::startOver() {
    mpzIndex = 0;
    dblIndex = 0;
    SetStartZ(myReps, freqs, z, IsComb, n, m, dblIndex,
              mpzIndex, IsRep, IsMult, IsGmp);
}

SEXP Combo::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0) &&
        CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return BasicVecReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextIter(freqs, z, n1, m1);
        return BasicVecReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP Combo::prevComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        return BasicVecReturn();
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 1)) {
        decrement(IsGmp, mpzIndex, dblIndex);
        prevIter(freqs, z, n1, m1);
        return BasicVecReturn();
    } else if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    } else {
        return R_NilValue;
    }
}

SEXP Combo::nextNumCombs(SEXP RNum) {

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
            nextIter(freqs, z, n1, m1);
        }

        return MatForward(nRows, numIncrement);
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP Combo::prevNumCombs(SEXP RNum) {

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
            prevIter(freqs, z, n1, m1);
        }

        decrement(IsGmp, mpzIndex, dblIndex, numDecrement);
        return MatReverse(nRows);
    } else if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    } else {
        return R_NilValue;
    }
}

SEXP Combo::nextGather() {

    if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        return ToSeeLast();
    }

    if (IsGmp) {
        mpzTemp = computedRowsMpz - mpzIndex;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = computedRows - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextIter(freqs, z, n1, m1);
        }

        if (IsGmp) {
            mpzIndex = computedRowsMpz + 1;
        } else {
            dblIndex = computedRows + 1;
        }

        return MatForward(nRows, 0);
    } else {
        return R_NilValue;
    }
}

SEXP Combo::prevGather() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 1)) {
        return ToSeeFirst();
    }

    if (IsGmp) {
        mpzTemp = mpzIndex - 1;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = dblIndex - 1;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                       computedRowsMpz, computedRows, true)) {
            prevIter(freqs, z, n1, m1);
        }

        if (IsGmp) {
            mpzIndex = 0;
        } else {
            dblIndex = 0;
        }

        return MatReverse(nRows);
    } else {
        return R_NilValue;
    }
}

SEXP Combo::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return BasicVecReturn();
    } else {
        return ToSeeFirst(false);
    }
}

SEXP Combo::randomAccess(SEXP RindexVec) {

    std::size_t sampSize;
    std::vector<double> mySample;
    SetIndexVec(RindexVec, mySample, sampSize, IsGmp, computedRows);

    const std::size_t bigSampSize = IsGmp ? sampSize : 1;
    std::vector<mpz_class> mpzVec(bigSampSize);

    if (IsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec, sampSize, computedRowsMpz);
    }

    if (sampSize > 1) {
        int nThreads = 1;
        bool LocalPar = Parallel;
        const int limit = 2;

        SetThreads(LocalPar, maxThreads, sampSize,
                   myType, nThreads, sexpNThreads, limit);

        return SampCombPermMain(sexpVec, vInt, vNum, mySample, mpzVec,
                                myReps, nthResFun, myType, n, m, sampSize,
                                nThreads, false, IsGmp, Parallel);
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
        return BasicVecReturn();
    }
}

SEXP Combo::front() {

    if (IsGmp) {
        mpzIndex = 1;
        mpzTemp  = 0;
    } else {
        dblIndex = 1;
        dblTemp  = 0;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return BasicVecReturn();
}

SEXP Combo::back() {

    if (IsGmp) {
        mpzIndex = computedRowsMpz;
        mpzTemp  = computedRowsMpz - 1;
    } else {
        dblIndex = computedRows;
        dblTemp  = computedRows - 1;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return BasicVecReturn();
}

SEXP Combo::sourceVector() const {
    return sexpVec;
}

SEXP Combo::summary() {
    const std::string CoPerm = IsComb ? "Combinations " : "Permutations ";
    const std::string RepStr = IsRep ? "with repetition " : "";
    const std::string MultiStr = IsMult ? "of a multiset " : "";
    const std::string strDesc = CoPerm + RepStr + MultiStr + "of "
                                + std::to_string(n) + " choose " +
                                    std::to_string(m);
    const double dblDiff = IsGmp ? 0 : computedRows - dblIndex;

    if (IsGmp) {
        mpzTemp = computedRowsMpz - mpzIndex;
    }

    const char *names[] = {"description", "currentIndex",
                           "totalResults", "totalRemaining", ""};

    cpp11::sexp res = Rf_mkNamed(VECSXP, names);

    SET_VECTOR_ELT(res, 0, Rf_mkString(strDesc.c_str()));
    SET_VECTOR_ELT(res, 1, CppConvert::GetCount(IsGmp, mpzIndex, dblIndex));
    SET_VECTOR_ELT(res, 2, CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows));
    SET_VECTOR_ELT(res, 3, CppConvert::GetCount(IsGmp, mpzTemp, dblDiff));
    return res;
}
