#include "ComboGroups/ComboGroupsClass.h"

SEXP ComboGroupsClass::SingleReturn() {

    cpp11::sexp res = BasicVecReturn();

    if (IsArray) {
        Rf_setAttrib(res, R_DimSymbol, dim);
        Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
    } else {
        Rf_setAttrib(res, R_NamesSymbol, myNames);
    }

    return res;
}

SEXP ComboGroupsClass::GeneralReturn(int numResults) {

    int nThreads = 1;
    bool LocalPar = Parallel;
    const int limit = 20000;
    std::vector<double> tempSample;
    std::vector<mpz_class> tempBigSamp;

    SetThreads(LocalPar, maxThreads, numResults,
               myType, nThreads, sexpNThreads, limit);

    cpp11::sexp res = GetComboGroups(
        sexpVec, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp, FinalTouch, vNum,
        vInt, z, myType, tempSample, tempBigSamp, mpzIndex, dblIndex, n,
        numResults, nThreads, IsArray, false, LocalPar, false, IsGmp
    );

    zUpdateIndex(vNum, vInt, z, sexpVec, res, m, numResults);
    return res;
}

ComboGroupsClass::ComboGroupsClass(
    SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
    const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
    const std::vector<int> &RvInt, const std::vector<double> &RvNum,
    VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
    SEXP RNumGroups, SEXP RGrpSize, SEXP RRetType
) : Combo(Rv, Rm, RcompRows, bVec, Rreps, Rfreqs, RvInt, RvNum, typePass,
          RmaxThreads, RnumThreads, Rparallel),
          CmbGrp(GroupPrep(Rv, RNumGroups, RGrpSize, n)) {

    prevIterAvailable = false;
    CmbGrp->SetCount();
    CmbGrpClsFuncs f = GetClassFuncs(CmbGrp);

    nthCmbGrp = f.nthDbl;
    nthCmbGrpGmp = f.nthGmp;
    nextCmbGrp = f.next;
    FinalTouch = f.finishing;

    IsGmp = CmbGrp->GetIsGmp();
    computedRows = CmbGrp->GetDblCount();
    if (IsGmp) computedRowsMpz = CmbGrp->GetMpzCount();

    z.resize(n);
    std::iota(z.begin(), z.end(), 0);

    std::string retType(CHAR(STRING_ELT(RRetType, 0)));

    if (retType != "3Darray" && retType != "matrix") {
        cpp11::stop("retType must be '3Darray' or 'matrix'");
    }

    if (retType == "3Darray" && CmbGrp->GetType() != "Uniform") {
        std::string msg = "3Darray output is not possible! Using matrix instead.";
        cpp11::message(msg.c_str());
        retType = "matrix";
    }

    IsArray = (retType == "3Darray");
    r = CmbGrp->GetNumGrps();
    rDisp = r;
    const int grpSize = n / r;

    std::vector<std::string> myColNames(r, "Grp");

    for (int j = 0; j < r; ++j) {
        myColNames[j] += std::to_string(j + 1);
    }

    for (auto g: CmbGrp->GetGroupSizes()) {
        grpSizeDesc += (std::to_string(g) + ", ");
    }

    if (IsArray) {
        myNames.resize(r);

        for (int i = 0; i < r; ++i) {
            myNames[i] = myColNames[i].c_str();
        }

        dimNames.resize(2);
        dimNames[1] = myNames;

        cpp11::integers temp_dim({grpSize, r});
        dim = temp_dim;
    } else if (CmbGrp->GetType() == "Uniform") {
        myNames.resize(n);

        for (int i = 0, k = 0; i < r; ++i) {
            for (int j = 0; j < grpSize; ++j, ++k) {
                myNames[k] = myColNames[i].c_str();
            }
        }
    } else if (CmbGrp->GetOneGrp()) {
        myNames.resize(n);
        std::vector<int> vGrpSizes(CmbGrp->GetGroupSizes());

        const int numOneGrps = vGrpSizes.front();
        std::vector<int> realGrps(vGrpSizes);
        realGrps.erase(realGrps.begin());
        realGrps.insert(realGrps.begin(), numOneGrps, 1);

        rDisp = realGrps.size();
        std::vector<std::string> myColNamesOne(rDisp, "Grp");

        for (int j = 0; j < rDisp; ++j) {
            myColNamesOne[j] += std::to_string(j + 1);
        }

        for (int i = 0, k = 0; i < rDisp; ++i) {
            for (int j = 0; j < realGrps[i]; ++j, ++k) {
                myNames[k] = myColNamesOne[i].c_str();
            }
        }

        grpSizeDesc.clear();

        for (auto g: realGrps) {
            grpSizeDesc += (std::to_string(g) + ", ");
        }
    } else {
        myNames.resize(n);
        std::vector<int> vGrpSizes(CmbGrp->GetGroupSizes());

        for (int i = 0, k = 0; i < r; ++i) {
            for (int j = 0; j < vGrpSizes[i]; ++j, ++k) {
                myNames[k] = myColNames[i].c_str();
            }
        }
    }

    // Remove the last space and comma
    grpSizeDesc.pop_back();
    grpSizeDesc.pop_back();
}

void ComboGroupsClass::startOver() {
    mpzIndex = 0;
    dblIndex = 0;
    std::iota(z.begin(), z.end(), 0);
}

SEXP ComboGroupsClass::nextComb() {

    if (CheckEqSi(IsGmp, mpzIndex, dblIndex, 0) &&
        CheckIndLT(IsGmp, mpzIndex, dblIndex,
                   computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        return SingleReturn();
    } else if (CheckIndLT(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        increment(IsGmp, mpzIndex, dblIndex);
        nextCmbGrp(z);
        return SingleReturn();
    } else if (CheckEqInd(IsGmp, mpzIndex, dblIndex,
                          computedRowsMpz, computedRows)) {
        return ToSeeLast();
    } else {
        return R_NilValue;
    }
}

SEXP ComboGroupsClass::nextNumCombs(SEXP RNum) {

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
            nextCmbGrp(z);
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

SEXP ComboGroupsClass::nextGather() {

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
            nextCmbGrp(z);
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

SEXP ComboGroupsClass::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex, dblIndex,
                    computedRowsMpz, computedRows)) {
        return ToSeeLast(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex, dblIndex, 0)) {
        return SingleReturn();
    } else {
        return ToSeeFirst(false);
    }
}

SEXP ComboGroupsClass::randomAccess(SEXP RindexVec) {

    std::size_t sampSize;
    std::vector<double> mySample;
    const bool SampIsGmp = (computedRows > SampleLimit);
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

        cpp11::sexp res = GetComboGroups(
            sexpVec, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp, FinalTouch, vNum,
            vInt, z, myType, mySample, mpzVec, mpzIndex, dblIndex, n,
            sampSize, nThreads, IsArray, false, LocalPar, true, IsGmp
        );

        z = before;
        return res;
    } else {
        if (IsGmp) {
            mpzIndex = mpzVec.front() + 1;
            mpzTemp  = mpzVec.front();
            z = nthCmbGrpGmp(mpzTemp);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp  = mySample.front();
            z = nthCmbGrp(dblTemp);
        }

        return SingleReturn();
    }
}

SEXP ComboGroupsClass::front() {

    if (IsGmp) {
        mpzIndex = 1;
        mpzTemp  = 0;
        z = nthCmbGrpGmp(mpzTemp);
    } else {
        dblIndex = 1;
        dblTemp  = 0;
        z = nthCmbGrp(dblTemp);
    }

    return SingleReturn();
}

SEXP ComboGroupsClass::back() {

    if (IsGmp) {
        mpzIndex = computedRowsMpz;
        mpzTemp  = computedRowsMpz - 1;
        z = nthCmbGrpGmp(mpzTemp);
    } else {
        dblIndex = computedRows;
        dblTemp  = computedRows - 1;
        z = nthCmbGrp(dblTemp);
    }

    return SingleReturn();
}

SEXP ComboGroupsClass::summary() {

    const std::string gtype = CmbGrp->GetType();
    const std::string prefix = "Partition of v of length " +
        std::to_string(n) + " into " + std::to_string(rDisp);
    const std::string suffix = (gtype == "Uniform") ? " uniform groups" :
        " groups of sizes: " + grpSizeDesc;

    const std::string strDesc = prefix + suffix;
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
