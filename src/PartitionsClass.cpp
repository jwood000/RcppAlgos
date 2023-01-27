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

    if (paragon) {
        for (auto &z_i: z) {
            z_i = vInt[z_i];
        }
    }

    SetPartValues();
}

SEXP Partitions::MultisetMatrix(int nRows) {

    cpp11::sexp res = Rf_allocMatrix(RTYPE, nRows, nCols);
    const std::size_t lastRow = nRows - 1;
    const std::size_t unRows  = nRows;
    const std::size_t unCols  = nCols;

    if (RTYPE == INTSXP) {
        int* ptrOut = INTEGER(res);

        for (std::size_t i = 0; i < lastRow; ++i) {
            for (std::size_t j = 0; j < unCols; ++j) {
                ptrOut[i + j * unRows] = vInt[z[j]];
            }

            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        for (std::size_t j = 0; j < unCols; ++j) {
            ptrOut[lastRow + j * unRows] = vInt[z[j]];
        }
    } else {
        double* ptrOut = REAL(res);

        for (std::size_t i = 0; i < lastRow; ++i) {
            for (std::size_t j = 0; j < unCols; ++j) {
                ptrOut[i + j * unRows] = vNum[z[j]];
            }

            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        for (std::size_t j = 0; j < unCols; ++j) {
            ptrOut[lastRow + j * unRows] = vNum[z[j]];
        }
    }

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
    bool RnumUnknown, double RcnstrtRows, const mpz_class &RcnstrtRowsMpz
) : ComboRes(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
             RmaxThreads, RnumThreads, Rparallel, Rpart, RcompVec, RtarVals,
             RtarIntVals, RstartZ, RmainFun, RFunTest, RfunDbl, Rctype,
             RstrtLen, Rcap, RKeepRes, RnumUnknown, RcnstrtRows,
             RcnstrtRowsMpz),
    paragon(ctype == ConstraintType::PartStandard),
    stdPartNext(paragon && !part.isComp),
    stdCompZeroSpesh(paragon && part.isComp && !part.isWeak),
    genCompZeroSpesh(!paragon && part.isComp && !part.isWeak &&
        part.includeZero),
    lastCol(part.width - 1), lastElem(n - 1),
    nextParts(
        GetNextPartsPtr(
            part.ptype,
            !(stdPartNext || stdCompZeroSpesh || genCompZeroSpesh),
            part.isComp
        )
    ),
    nthParts((part.ptype == PartitionType::LengthOne ||
              part.ptype == PartitionType::Multiset  ||
              CheckEqSi(part.isGmp, cnstrtCountMpz, cnstrtCount, 0)) ?
              nullptr : GetNthPartsFunc(part.ptype, part.isGmp,
                                        part.isComp)) {

    bAddOne = paragon && !part.includeZero;
    rpsCnt = myReps;
    IsGmp  = part.isGmp;
    SetPartValues();
    prevIterAvailable = false;
}

void Partitions::startOver() {
    mpzIndex = 0;
    dblIndex = 0;
    rpsCnt = myReps;
    z = part.startZ;
    SetPartValues();
}

SEXP Partitions::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0) &&
        CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {
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
    CppConvert::convertPrimitive(RNum, num, VecType::Integer,
                                   "The number of results");

    if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   cnstrtCountMpz, cnstrtCount)) {

        int nRows = 0;
        int numIncrement = 0;

        if (IsGmp) {
            mpzTemp = cnstrtCountMpz - mpzIndex;
            nRows = cmp(mpzTemp, num) < 0 ? mpzTemp.get_si() : num;
            numIncrement = cmp(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
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
            bUpper = true;
            cpp11::sexp res = MatrixReturn(nRows);
            increment(IsGmp, mpzIndex, dblIndex, numIncrement);
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows, bAddOne);
            SetPartValues();
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
        mpzTemp = cnstrtCountMpz - mpzIndex;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
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

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextParts(rpsCnt, z, edge, boundary, pivot,
                      tarDiff, lastCol, lastElem);
        }

        if (IsGmp) {
            mpzIndex = cnstrtCountMpz + 1;
        } else {
            dblIndex = cnstrtCount + 1;
        }

        if (part.ptype == PartitionType::Multiset) {
            return MultisetMatrix(nRows);
        } else {
            bUpper = false;
            cpp11::sexp res = MatrixReturn(nRows);
            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, nRows, bAddOne);
            SetPartValues();
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

    const std::size_t bigSampSize = SampIsGmp ? sampSize : 1;
    std::vector<mpz_class> mpzVec(bigSampSize);

    if (SampIsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec, sampSize, cnstrtCountMpz);
    }

    if (sampSize > 1) {
        int nThreads = 1;
        bool LocalPar = Parallel;
        const int limit = 2;

        SetThreads(LocalPar, maxThreads, sampSize,
                   myType, nThreads, sexpNThreads, limit);

        if (myType == VecType::Integer) {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, sampSize, part.width);
            int* matInt = INTEGER(res);

            ThreadSafeSample(matInt, res, vInt, mySample, mpzVec,
                             myReps, nthParts, part.width, sampSize,
                             nThreads, Parallel, false, part.mapTar,
                             strtLen, cap, IsGmp);

            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, sampSize, bAddOne);
            SetPartValues();
            return res;
        } else {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, sampSize, part.width);
            double* matNum = REAL(res);

            ThreadSafeSample(matNum, res, vNum, mySample, mpzVec,
                             myReps, nthParts, part.width, sampSize,
                             nThreads, Parallel, false, part.mapTar,
                             strtLen, cap, IsGmp);

            zUpdateIndex(vNum, vInt, z, sexpVec, res, width, sampSize, bAddOne);
            SetPartValues();
            return res;
        }
    } else {
        if (IsGmp) {
            mpzIndex = mpzVec.front() + 1;
            mpzTemp  = mpzVec.front();
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp  = mySample.front();
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
        mpzIndex = 1;
        mpzTemp  = 0;
    } else {
        dblIndex = 1;
        dblTemp  = 0;
    }

    MoveZToIndex();
    return VecReturn();
}

SEXP Partitions::back() {

    if (nthParts == nullptr) {
        cpp11::stop("No random access available for this scenario");
    }

    if (IsGmp) {
        mpzIndex = cnstrtCountMpz;
        mpzTemp  = cnstrtCountMpz - 1;
    } else {
        dblIndex = cnstrtCount;
        dblTemp  = cnstrtCount - 1;
    }

    MoveZToIndex();
    return VecReturn();
}

SEXP Partitions::summary() {
    const std::string RepStr = (IsRep) ? "with repetition " : "";
    const std::string MultiStr = (IsMult) ? "of a multiset " : "";
    const std::string strDesc = (part.isComp ? "Compositions " : "Partitions ")
          + RepStr + MultiStr + "of " + std::to_string(part.target) +
              " into " + std::to_string(width) + " parts";
    const double dblDiff = IsGmp ? 0 : cnstrtCount - dblIndex;

    if (IsGmp) mpzTemp = cnstrtCountMpz - mpzIndex;
    const char *names[] = {"description", "currentIndex",
                           "totalResults", "totalRemaining", ""};

    cpp11::sexp res = Rf_mkNamed(VECSXP, names);

    SET_VECTOR_ELT(res, 0, Rf_mkString(strDesc.c_str()));
    SET_VECTOR_ELT(res, 1, CppConvert::GetCount(IsGmp, mpzIndex, dblIndex));
    SET_VECTOR_ELT(res, 2, CppConvert::GetCount(IsGmp, cnstrtCountMpz, cnstrtCount));
    SET_VECTOR_ELT(res, 3, CppConvert::GetCount(IsGmp, mpzTemp, dblDiff));
    return res;
}
