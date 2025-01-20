#include "Cartesian/CartesianClass.h"

// These 3 previous iterators methods are here because we are because we
// are inheriting from the parent Iterator class, the class that is utilized
// in ExposeClass.cpp to connect the S4 class methods on the R side to C++.
// ****************************************************************************
SEXP CartesianClass::prevIter() {
    cpp11::stop("No prevIter available yet for Cartesian Class");
}

SEXP CartesianClass::prevNumIters(SEXP RNum) {
    cpp11::stop("No prevNumIters available yet for Cartesian Class");
}

SEXP CartesianClass::prevGather() {
    cpp11::stop("No prevGather available yet for Cartesian Class");
}
// ****************************************************************************

SEXP CartesianClass::VectorReturn() {
    switch (myType) {
        case VecType::Logical : {
            cpp11::sexp res = Rf_allocVector(LGLSXP, nCols);
            int* ptrOut = LOGICAL(res);

            for (int j = 0; j < nCols; ++j) {
                ptrOut[j] = boolVec[idx[j + z[j]]];
            }

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocVector(INTSXP, nCols);
            int* ptrOut = INTEGER(res);

            for (int j = 0; j < nCols; ++j) {
                ptrOut[j] = intVec[idx[j + z[j]]];
            }

            if (typeCheck[tFac]) SetFactorClass(res, RList[0]);
            return res;
        } case VecType::Character : {
            cpp11::sexp res = Rf_allocVector(STRSXP, nCols);

            for (int j = 0; j < nCols; ++j) {
                SET_STRING_ELT(res, j, STRING_ELT(charVec, idx[j + z[j]]));
            }

            return res;
        } case VecType::Complex : {
            cpp11::sexp res = Rf_allocVector(CPLXSXP, nCols);
            Rcomplex* ptrOut = COMPLEX(res);

            for (int j = 0; j < nCols; ++j) {
                ptrOut[j] = cmplxVec[idx[j + z[j]]];
            }

            return res;
        } case VecType::Raw : {
            cpp11::sexp res = Rf_allocVector(RAWSXP, nCols);
            Rbyte* ptrOut = RAW(res);

            for (int j = 0; j < nCols; ++j) {
                ptrOut[j] = rawVec[idx[j + z[j]]];
            }

            return res;
        } case VecType::Numeric : {
            cpp11::sexp res = Rf_allocVector(REALSXP, nCols);
            double* ptrOut = REAL(res);

            for (int j = 0; j < nCols; ++j) {
                ptrOut[j] = dblVec[idx[j + z[j]]];
            }

            return res;
        } default : {
            cpp11::stop("Only atomic types are supported for v");
        }
    }
}

SEXP CartesianClass::SingleReturn() {

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);

        for (int j = 0; j < nCols; ++j) {
            switch (TYPEOF(RList[j])) {
                case INTSXP: {
                    cpp11::sexp res = Rf_allocVector(INTSXP, 1);
                    int* intSexpVec = INTEGER(res);
                    intSexpVec[0] = intVec[idx[j + z[j]]];
                    if (IsFactor[j]) SetFactorClass(res, RList[j]);
                    DataFrame[j] = res;
                    break;
                } case LGLSXP: {
                    cpp11::sexp res = Rf_allocVector(LGLSXP, 1);
                    int* boolSexpVec = LOGICAL(res);
                    boolSexpVec[0] = boolVec[idx[j + z[j]]];
                    DataFrame[j] = res;
                    break;
                } case REALSXP: {
                    cpp11::sexp res = Rf_allocVector(REALSXP, 1);
                    double* dblSexpVec = REAL(res);
                    dblSexpVec[0] = dblVec[idx[j + z[j]]];
                    DataFrame[j] = res;
                    break;
                } case CPLXSXP : {
                    cpp11::sexp res = Rf_allocVector(CPLXSXP, 1);
                    Rcomplex* cmplxSexpVec = COMPLEX(res);
                    cmplxSexpVec[0] = cmplxVec[idx[j + z[j]]];
                    DataFrame[j] = res;
                    break;
                } case RAWSXP : {
                    cpp11::sexp res = Rf_allocVector(RAWSXP, 1);
                    Rbyte* rawSexpVec = RAW(res);
                    rawSexpVec[0] = rawVec[idx[j + z[j]]];
                    DataFrame[j] = res;
                    break;
                } case STRSXP: {
                    cpp11::sexp res = Rf_allocVector(STRSXP, 1);
                    SET_STRING_ELT(res, 0, STRING_ELT(charVec, idx[j + z[j]]));
                    DataFrame[j] = res;
                    break;
                } default : {
                    cpp11::stop("Only atomic types are supported for v");
                }
            }
        }

        DataFrame.attr("row.names") = {1};
        DataFrame.names() = RList.names();
        DataFrame.attr("class") = "data.frame";
        return DataFrame;
    } else {
        cpp11::sexp res = VectorReturn();
        res.names() = RList.names();
        return res;
    }
}

SEXP CartesianClass::GeneralReturn(int numResults) {

    int nThreads = 1;
    bool LocalPar = Parallel;
    const int limit = 20000;
    std::vector<double> tempSample;
    std::vector<mpz_class> tempBigSamp;

    SetThreads(LocalPar, maxThreads, numResults,
               myType, nThreads, sexpNThreads, limit);

    cpp11::sexp res = GetProduct(
        idx, typeCheck, IsFactor, RList, intVec, dblVec,
        boolVec, cmplxVec, rawVec, charVec, lenGrps, z,
        tempSample, tempBigSamp, dblIndex, mpzIndex, numResults,
        nCols, IsDF, nThreads, LocalPar, IsGmp, false
    );

    mpzTemp = mpzIndex - 1;
    dblTemp = dblIndex - 1;
    GetStartProd(lenNxtPr, z, mpzTemp, dblTemp, 0, IsGmp);

    SetMatrixColnames(res, RList.names());
    return res;
}

CartesianClass::CartesianClass(
    SEXP Rv_RList, SEXP RcompRows, int RmaxThreads, SEXP RnumThreads,
    bool Rparallel, bool RIsGmp, const std::vector<int> &Ridx,
    const std::vector<int> &RtypeCheck, const std::vector<int> &RIsFactor,
    const std::vector<int> &RintVec, const std::vector<double> &RdblVec,
    const std::vector<int> &RlglVec, const std::vector<Rcomplex> &RcplxVec,
    const std::vector<Rbyte> &RrawVec, const cpp11::strings &RcharVec,
    const std::vector<int> &RlenGrps, bool RisDF, int RnCols, VecType RmyType
) : Iterator(Rv_RList, VecType::Raw, RcompRows, RmaxThreads,
             RnumThreads, Rparallel, RIsGmp), RList(Rv_RList), idx(Ridx),
    lenGrps(RlenGrps), typeCheck(RtypeCheck), IsFactor(RIsFactor),
    intVec(RintVec), dblVec(RdblVec), boolVec(RlglVec), cmplxVec(RcplxVec),
    rawVec(RrawVec), charVec(RcharVec), IsDF(RisDF),
    nCols(RnCols), myType(RmyType) {

    dblIndex = 0;
    mpzIndex = 0;

    z.assign(nCols, 0);
    lenNxtPr = lenGrps;

    for (auto &v_i: lenNxtPr) {
        v_i = (v_i / nCols) + 1;
    }
}

void CartesianClass::startOver() {
    mpzIndex = 0;
    dblIndex = 0;
    std::fill(z.begin(), z.end(), 0);
}

SEXP CartesianClass::nextIter() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0) &&
        CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return SingleReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextProduct(lenGrps, z, nCols);
        return SingleReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP CartesianClass::nextNumIters(SEXP RNum) {

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
            nextProduct(lenGrps, z, nCols);
        }

        increment(IsGmp, mpzIndex, dblIndex, numIncrement);
        return GeneralReturn(nRows);
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP CartesianClass::nextGather() {

    if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        return ToSeeLast();
    }

    if (IsGmp) {
        mpzTemp = computedRowsMpz - mpzIndex;

        if (cmp(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            cpp11::stop(
                "The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str()
            );
        }
    } else {
        dblTemp = computedRows - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            cpp11::stop(
                "The number of requested rows is greater than %s",
                std::to_string(std::numeric_limits<int>::max()).c_str()
            );
        }
    }

    const int nRows = IsGmp ? mpzTemp.get_si() : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
            nextProduct(lenGrps, z, nCols);
        }

        if (IsGmp) {
            mpzIndex = computedRowsMpz + 1;
        } else {
            dblIndex = computedRows + 1;
        }

        return GeneralReturn(nRows);
    } else {
        return R_NilValue;
    }
}

SEXP CartesianClass::currIter() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return SingleReturn();
    } else {
        return ToSeeFirst(false);
    }
}

SEXP CartesianClass::randomAccess(SEXP RindexVec) {

    std::size_t sampSize;
    std::vector<double> mySample;
    const bool SampIsGmp = IsGmp || computedRows > SampleLimit;
    SetIndexVec(RindexVec, mySample, sampSize, SampIsGmp, computedRows);

    const std::size_t bigSampSize = SampIsGmp ? sampSize : 1;
    std::vector<mpz_class> mpzVec(bigSampSize);

    if (SampIsGmp) {
        SetIndexVecMpz(RindexVec, mpzVec, sampSize, computedRowsMpz);
    }

    if (sampSize > 1) {
        int nThreads = 1;
        bool LocalPar = Parallel;
        const int limit = 2;

        SetThreads(LocalPar, maxThreads, sampSize,
                   myType, nThreads, sexpNThreads, limit);

        const std::vector<int> before(z);

        cpp11::sexp res = GetProduct(
            idx, typeCheck, IsFactor, RList, intVec, dblVec,
            boolVec, cmplxVec, rawVec, charVec, lenGrps, z,
            mySample, mpzVec, dblIndex, mpzIndex, sampSize,
            nCols, IsDF, nThreads, LocalPar, IsGmp, true
        );

        z = before;

        SetMatrixColnames(res, RList.names());
        return res;
    } else {
        if (IsGmp) {
            mpzIndex = mpzVec.front() + 1;
            mpzTemp  = mpzVec.front();
            z = nthProductGmp(mpzTemp, lenNxtPr);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp  = mySample.front();
            z = nthProduct(dblTemp, lenNxtPr);
        }

        return SingleReturn();
    }
}

SEXP CartesianClass::front() {

    if (IsGmp) {
        mpzIndex = 1;
        mpzTemp  = 0;
    } else {
        dblIndex = 1;
        dblTemp  = 0;
    }

    std::fill(z.begin(), z.end(), 0);
    return SingleReturn();
}

SEXP CartesianClass::back() {

    if (IsGmp) {
        mpzIndex = computedRowsMpz;
        mpzTemp  = computedRowsMpz - 1;
    } else {
        dblIndex = computedRows;
        dblTemp  = computedRows - 1;
    }

    GetStartProd(lenNxtPr, z, mpzTemp, dblTemp, 0, IsGmp);
    return SingleReturn();
}

SEXP CartesianClass::summary() {

    std::string basic = "Cartesian Product of the source";
    std::string mInfo = "see the sourceVector method for more info";
    const std::string strDesc = basic + " (" + mInfo + ")";
    const double dblDiff = IsGmp ? 0 : computedRows - dblIndex;

    if (IsGmp) mpzTemp = computedRowsMpz - mpzIndex;
    const char *names[] = {"description", "currentIndex",
                           "totalResults", "totalRemaining", ""};

    cpp11::sexp res = Rf_mkNamed(VECSXP, names);

    SET_VECTOR_ELT(res, 0, Rf_mkString(strDesc.c_str()));
    SET_VECTOR_ELT(res, 1, CppConvert::GetCount(IsGmp, mpzIndex, dblIndex));
    SET_VECTOR_ELT(res, 2, CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows));
    SET_VECTOR_ELT(res, 3, CppConvert::GetCount(IsGmp, mpzTemp, dblDiff));
    return res;
}
