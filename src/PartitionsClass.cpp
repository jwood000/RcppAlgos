#include "Partitions/PartitionsClass.h"

void Partitions::SetPartValues() {
    if (part.ptype == PartitionType::Multiset) {
        PrepareMultisetPart(rpsCnt, z, boundary, pivot,
                            edge, lastCol, lastElem);
    } else if (std::find(RepPTypeArr.cbegin(), RepPTypeArr.cend(),
                         part.ptype) != RepPTypeArr.cend()) {
        PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);
    } else {
        PrepareDistinctPart(z, boundary, pivot, edge,
                            tarDiff, lastElem, lastCol);
    }
}

void Partitions::MoveZToIndex() {
    z = nthParts(part.mapTar, width, cap, strtLen, dblTemp, mpzTemp);

    if (ctype == ConstraintType::PartStandard) {
        for (auto &z_i: z) {
            z_i = vInt[z_i];
        }
    }

    SetPartValues();
}

SEXP Partitions::MultisetMatrix(int nRows) {

    SEXP res = PROTECT(Rf_allocMatrix(RTYPE, nRows, nCols));
    const int lastRow = nRows - 1;

    if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);

        for (int i = 0; i < lastRow; ++i) {
            for (int j = 0; j < nCols; ++j) {
                ptrOut[i + j * nRows] = vInt[z[j]];
            }

            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        for (int j = 0; j < nCols; ++j) {
            ptrOut[lastRow + j * nRows] = vInt[z[j]];
        }
    } else {
        double* ptrOut = REAL(res);

        for (int i = 0; i < lastRow; ++i) {
            for (int j = 0; j < nCols; ++j) {
                ptrOut[i + j * nRows] = vNum[z[j]];
            }

            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        for (int j = 0; j < nCols; ++j) {
            ptrOut[lastRow + j * nRows] = vNum[z[j]];
        }
    }

    UNPROTECT(1);
    return res;
}

Partitions::Partitions(
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
) : ComboRes(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
             RmaxThreads, RnumThreads, Rparallel, Rpart, RcompVec, RtarVals,
             RtarIntVals, RstartZ, RmainFun, RFunTest, RfunDbl, Rctype,
             RstrtLen, Rcap, RKeepRes, RnumUnknown, RcnstrtRows,
             RcnstrtRowsMpz),
    lastCol(part.width - 1), lastElem(n - 1),
    nextParts(GetNextPartsPtr(part.ptype,
                              ctype != ConstraintType::PartStandard)),
    nthParts(part.ptype == PartitionType::Multiset ? nullptr :
               GetNthPartsFunc(part.ptype, part.isGmp)) {

    bAddOne = (ctype == ConstraintType::PartStandard) && !part.includeZero;
    rpsCnt = myReps;
    IsGmp  = part.isGmp;
    SetPartValues();
    prevIterAvailable = false;
}

void Partitions::startOver() {
    if (IsGmp) {
      mpz_set_ui(mpzIndex, 0u);
    } else {
      dblIndex = 0;
    }

    rpsCnt = myReps;
    z = part.startZ;
    SetPartValues();
}

SEXP Partitions::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return VecReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextParts(rpsCnt, z, edge, boundary, pivot,
                  tarDiff, lastCol, lastElem);
        return VecReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP Partitions::nextNumCombs(SEXP RNum) {

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
            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        if (part.ptype == PartitionType::Multiset) {
            increment(IsGmp, mpzIndex, dblIndex, numIncrement);
            return MultisetMatrix(nRows);
        } else {
            bUpper   = true;
            SEXP res = PROTECT(MatrixReturn(nRows));
            increment(IsGmp, mpzIndex, dblIndex, numIncrement);
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows, bAddOne);
            SetPartValues();
            UNPROTECT(1);
            return res;
        }
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP Partitions::nextGather() {

    if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast();
    }

    if (IsGmp) {
        mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);

        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    } else {
        dblTemp = cnstrtCount - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str());
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        if (IsGmp) {
            mpz_add_ui(mpzIndex, cnstrtCountMpz, 1u);
        } else {
            dblIndex = cnstrtCount + 1;
        }

        if (part.ptype == PartitionType::Multiset) {
            return MultisetMatrix(nRows);
        } else {
            bUpper   = false;
            SEXP res = PROTECT(MatrixReturn(nRows));
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows, bAddOne);
            SetPartValues();
            UNPROTECT(1);
            return res;
        }
    } else {
        return R_NilValue;
    }
}

SEXP Partitions::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    cnstrtCountMpz, cnstrtCount)) {
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return VecReturn();
    } else {
        return ToSeeFirst(false);
    }
}

SEXP Partitions::randomAccess(SEXP RindexVec) {

    if (nthParts == nullptr) {
        cpp11::stop("No random access available for this scenario");
    }

    std::size_t sampSize;
    std::vector<double> mySample;
    const bool SampIsGmp = (cnstrtCount > SampleLimit);
    SetIndexVec(RindexVec, mySample, sampSize, SampIsGmp, cnstrtCount);

    const std::size_t bigSampSize = (SampIsGmp) ? sampSize : 1;
    auto mpzVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);

    for (std::size_t i = 0; i < bigSampSize; ++i) {
        mpz_init(mpzVec[i]);
    }

    if (SampIsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec.get(), sampSize, cnstrtCountMpz);
    }

    if (sampSize > 1) {
        int nThreads = 1;
        bool LocalPar = Parallel;
        const int limit = 2;

        SetThreads(LocalPar, maxThreads, sampSize,
                   myType, nThreads, sexpNThreads, limit);

        if (myType == VecType::Integer) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, sampSize, part.width));
            int* matInt = INTEGER(res);

            ThreadSafeSample(matInt, res, vInt, mySample, mpzVec.get(),
                             myReps, nthParts, part.width, sampSize,
                             nThreads, Parallel, false, part.mapTar,
                             strtLen, cap, IsGmp);
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, sampSize, bAddOne);
            SetPartValues();
            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, sampSize, part.width));
            double* matNum = REAL(res);

            ThreadSafeSample(matNum, res, vNum, mySample, mpzVec.get(),
                             myReps, nthParts, part.width, sampSize,
                             nThreads, Parallel, false, part.mapTar,
                             strtLen, cap, IsGmp);
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, sampSize, bAddOne);
            SetPartValues();
            UNPROTECT(1);
            return res;
        }
    } else {
        if (IsGmp) {
            mpz_add_ui(mpzIndex, mpzVec[0], 1u);
            mpz_set(mpzTemp, mpzVec[0]);
            mpz_clear(mpzVec[0]);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp = mySample.front();
        }

        MoveZToIndex();
        return VecReturn();
    }
}

SEXP Partitions::front() {

    if (nthParts == nullptr) {
        cpp11::stop("No random access available for this scenario");
    }

    if (IsGmp) {
        mpz_set_ui(mpzIndex, 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }

    MoveZToIndex();
    return VecReturn();
}

SEXP Partitions::back() {

    if (nthParts == nullptr) {
        cpp11::stop("No random access available for this scenario");
    }

    if (IsGmp) {
        mpz_set(mpzIndex, cnstrtCountMpz);
        mpz_sub_ui(mpzTemp, cnstrtCountMpz, 1u);
    } else {
        dblIndex = cnstrtCount;
        dblTemp = cnstrtCount - 1;
    }

    MoveZToIndex();
    return VecReturn();
}

SEXP Partitions::summary() {
    const std::string RepStr = (IsRep) ? "with repetition " : "";
    const std::string MultiStr = (IsMult) ? "of a multiset " : "";
    const std::string strDesc = "Partitions " + RepStr + MultiStr + "of "
          + std::to_string(part.target) + " into " + std::to_string(width) + " parts";
    const double dblDiff = (IsGmp) ? 0 : cnstrtCount - dblIndex;

    if (IsGmp) {
        mpz_sub(mpzTemp, cnstrtCountMpz, mpzIndex);
    }

    const char *names[] = {"description", "currentIndex",
                           "totalResults", "totalRemaining", ""};

    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));

    SET_VECTOR_ELT(res, 0, Rf_mkString(strDesc.c_str()));
    SET_VECTOR_ELT(res, 1, CleanConvert::GetCount(IsGmp, mpzIndex, dblIndex));
    SET_VECTOR_ELT(res, 2, CleanConvert::GetCount(IsGmp, cnstrtCountMpz, cnstrtCount));
    SET_VECTOR_ELT(res, 3, CleanConvert::GetCount(IsGmp, mpzTemp, dblDiff));

    UNPROTECT(1);
    return res;
}
