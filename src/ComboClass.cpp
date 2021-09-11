#include "ClassUtils/ComboClass.h"

SEXP Combo::VecReturn() {

    SEXP res = PROTECT(Rf_allocVector(RTYPE, m));

    switch (RTYPE) {
        case LGLSXP: {
            int* ptrOut = INTEGER(res);
            int* ptrIn  = INTEGER(sexpVec);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }

            break;
        } case INTSXP: {
            int* ptrOut = INTEGER(res);
            int* ptrIn  = INTEGER(sexpVec);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
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
            double* ptrIn  = REAL(sexpVec);

            for (int j = 0; j < m; ++j) {
                ptrOut[j] = ptrIn[z[j]];
            }

            break;
        }
    }

    UNPROTECT(1);
    return res;
}

SEXP Combo::MatForward(int nRows) {

    int nThreads = 1;
    bool LocalPar = Parallel;
    const int limit = 20000;
    
    SetThreads(LocalPar, maxThreads, nRows,
               myType, nThreads, SexpNThreads, limit);

    SEXP res = PROTECT(GetCombPerms(sexpVec, vNum, vInt, n, m, 0, true,
                                    IsComb, LocalPar, IsRep, IsMult, IsGmp,
                                    freqs, z, myReps, dblIndex, mpzIndex[0],
                                    nRows, nThreads, myType));

    zUpdateIndex(z, sexpVec, res, m, nRows);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    UNPROTECT(1);
    return res;
}

SEXP Combo::MatReverse(int nRows) {
    return GetPrevCombPerms(sexpVec, vNum, vInt, myReps, freqs, z,
                            prevIter, n, m, IsComb, IsMult, nRows, myType);
}

// The bVec Vector represents IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
Combo::Combo(SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
             const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
             const std::vector<int> &RvInt, const std::vector<double> &RvNum,
             VecType typePass, int RmaxThreads, SEXP RnumThreads)
            : SexpNThreads(RnumThreads), n(Rf_length(Rv)), m(Rm),
              m1(Rm - 1), RTYPE(TYPEOF(Rv)), maxThreads(RmaxThreads),
              sexpVec(Rv), IsFactor(bVec[0]), IsComb(bVec[1]),
              IsMult(bVec[2]), IsRep(bVec[3]), IsGmp(bVec[4]),
              Parallel(maxThreads > 1), computedRows(IsGmp ? 0 :
              Rf_asReal(RcompRows)), myType(typePass), vInt(RvInt),
              vNum(RvNum), freqs(Rfreqs), myReps(Rreps),
              n1(IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1)),
              myClass(bVec[0] ? Rf_getAttrib(Rv, R_ClassSymbol) :
                          Rf_allocVector(STRSXP, 0)),
              myLevels(bVec[0] ? Rf_getAttrib(Rv, R_LevelsSymbol) :
                           R_NilValue),
              nthResFun(GetNthResultFunc(bVec[1], bVec[2], bVec[3], bVec[4])),
              nextIter(GetNextIterPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
              prevIter(GetPrevIterPtr(bVec[1], bVec[2], bVec[3], bVec[5])) {

    z.resize(Rm);

    // Initialize trivial mpz_t value for functions which require mpz_t
    mpz_init(mpzTemp);
    mpz_init(mpzIndex[0]);
    mpz_init(computedRowsMpz[0]);

    if (IsGmp) {
        createMPZArray(RcompRows, computedRowsMpz, 1, "index vector");
    }

    dblIndex = 0;
    
    SetStartZ(myReps, freqs, z, IsComb, n, m, dblIndex,
              mpzIndex[0], IsRep, IsMult, IsGmp);
}

void Combo::startOver() {
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 0u);
    } else {
        dblIndex = 0;
    }
    
    SetStartZ(myReps, freqs, z, IsComb, n, m, dblIndex,
              mpzIndex[0], IsRep, IsMult, IsGmp);
}

SEXP Combo::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        return VecReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex,
                          computedRowsMpz[0], computedRows)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        nextIter(freqs, z, n1, m1);
        return VecReturn();
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

SEXP Combo::prevComb() {
    
    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return VecReturn();
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        prevIter(freqs, z, n1, m1);
        return VecReturn();
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

SEXP Combo::nextNumCombs(SEXP RNum) {

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
        return MatForward(nRows);
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

SEXP Combo::prevNumCombs(SEXP RNum) {

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

        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);

        decrement(IsGmp, mpzIndex[0], dblIndex, numDecrement);
        return MatReverse(nRows);
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

SEXP Combo::nextGather() {

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

        return MatForward(nRows);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP Combo::prevGather() {

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

        return MatReverse(nRows);
    } else {
        return Rf_ScalarLogical(false);
    }
}

SEXP Combo::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex,
                    computedRowsMpz[0], computedRows)) {
        const std::string message = "No more results. To see the last "
                                    "result, use the prevIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        return VecReturn();
    } else {
        const std::string message = "Iterator Initialized. To see the first "
                                    "result, use the nextIter method(s)\n\n";
        Rprintf(message.c_str());
        return Rf_ScalarLogical(0);
    }
}

SEXP Combo::randomAccess(SEXP RindexVec) {

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
        int nThreads = 1;
        bool LocalPar = Parallel;
        const int limit = 2;
        
        SetThreads(LocalPar, maxThreads, sampSize,
                   myType, nThreads, SexpNThreads, limit);

        return SampCombPermMain(sexpVec, vInt, vNum, mySample, mpzVec.get(),
                                myReps, nthResFun, myType, n, m, sampSize,
                                nThreads, false, IsGmp, Parallel);
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
        return VecReturn();
    }
}

SEXP Combo::front() {

    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecReturn();
}

SEXP Combo::back() {

    if (IsGmp) {
        mpz_set(mpzIndex[0], computedRowsMpz[0]);
        mpz_sub_ui(mpzTemp, computedRowsMpz[0], 1u);
    } else {
        dblIndex = computedRows;
        dblTemp = computedRows - 1;
    }

    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    if (!IsComb) TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    return VecReturn();
}

SEXP Combo::sourceVector() const {
    return sexpVec;
}

SEXP Combo::summary() {
    const std::string CoPerm = (IsComb) ? "Combinations " : "Permutations ";
    const std::string RepStr = (IsRep) ? "with repetition " : "";
    const std::string MultiStr = (IsMult) ? "of a multiset " : "";
    const std::string strDesc = CoPerm + RepStr + MultiStr + "of "
                                + std::to_string(n) + " choose " + std::to_string(m);
    const double dblDiff = (IsGmp) ? 0 : computedRows - dblIndex;

    if (IsGmp) {
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
    }

    const char *names[] = {"description", "currentIndex",
                           "totalResults", "totalRemaining", ""};

    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));

    SET_VECTOR_ELT(res, 0, Rf_mkString(strDesc.c_str()));
    SET_VECTOR_ELT(res, 1, CleanConvert::GetCount(IsGmp, mpzIndex[0], dblIndex));
    SET_VECTOR_ELT(res, 2, CleanConvert::GetCount(IsGmp, computedRowsMpz[0], computedRows));
    SET_VECTOR_ELT(res, 3, CleanConvert::GetCount(IsGmp, mpzTemp, dblDiff));

    UNPROTECT(1);
    return res;
}
