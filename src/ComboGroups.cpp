#include "ComboGroups/GetComboGroups.h"

[[cpp11::register]]
SEXP ComboGroupsCpp(SEXP Rv, SEXP RNumGroups, SEXP RGrpSize, SEXP RRetType,
                    SEXP Rlow, SEXP Rhigh, SEXP Rparallel, SEXP RNumThreads,
                    SEXP RmaxThreads, SEXP RIsSample, SEXP RindexVec,
                    SEXP RmySeed, SEXP RNumSamp, SEXP baseSample,
                    SEXP RNamed, SEXP myEnv) {

    int n;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    bool IsSample = CppConvert::convertFlag(RIsSample, "IsSample");
    bool Parallel = CppConvert::convertFlag(Rparallel, "Parallel");
    bool IsNamed  = (IsSample) ?
        CppConvert::convertFlag(RNamed, "namedSample") : false;

    std::vector<int> vInt;
    std::vector<double> vNum;

    SetType(myType, Rv);
    SetBasic(Rv, vNum, vInt, n, myType);

    const std::unique_ptr<ComboGroupsTemplate> CmbGrp =
        GroupPrep(Rv, RNumGroups, RGrpSize, n);

    CmbGrp->SetCount();
    const bool IsGmp = CmbGrp->GetIsGmp();

    if (myType == VecType::Integer) {
        vInt.assign(vNum.cbegin(), vNum.cend());
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    if (!IsSample) {
        SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper, lowerMpz,
                  upperMpz, CmbGrp->GetMpzCount(), CmbGrp->GetDblCount());
    }

    if (!IsGmp) lowerMpz = lower;
    std::vector<int> startZ;

    if (bLower && cmp(lowerMpz, 0) > 0) {
        startZ = IsGmp ? CmbGrp->nthComboGroupGmp(lowerMpz) :
                         CmbGrp->nthComboGroup(lower);
    } else {
        startZ.resize(n);
        std::iota(startZ.begin(), startZ.end(), 0);
    }

    int nRows = 0;

    if (!IsSample) {
        double userNumRows = 0;
        SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                      lowerMpz, lower, upper, CmbGrp->GetDblCount(),
                      CmbGrp->GetMpzCount(), nRows, userNumRows);
    }

    std::string retType(CHAR(STRING_ELT(RRetType, 0)));

    if (retType != "3Darray" && retType != "matrix") {
        cpp11::stop("retType must be '3Darray' or 'matrix'");
    }

    if (retType == "3Darray" && CmbGrp->GetType() != "Uniform") {
        std::string msg = "3Darray output is not possible! Using matrix instead.";
        cpp11::message(msg.c_str());
        retType = "matrix";
    }

    const bool IsArray = (retType == "3Darray");

    int sampSize;
    std::vector<double> mySample;

    if (IsSample) {
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                        CmbGrp->GetDblCount(), mySample, baseSample, myEnv);
    }

    const int bigSampSize = (IsSample && IsGmp) ? sampSize : 1;
    std::vector<mpz_class> myBigSamp(bigSampSize);

    if (IsSample) {
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                           IsGmp, CmbGrp->GetMpzCount(), myBigSamp);
    }

    const int numResults = (IsSample) ? sampSize : nRows;
    const int limit = (IsSample) ? 2 : 20000;
    SetThreads(Parallel, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    CmbGrpClsFuncs f = GetClassFuncs(CmbGrp);

    cpp11::sexp res = GetComboGroups(
        Rv, f.next, f.nthDbl, f.nthGmp, f.finishing, vNum, vInt,
        startZ, myType, mySample, myBigSamp, lowerMpz, lower, n, numResults,
        nThreads, IsArray, IsNamed, Parallel, IsSample, IsGmp
    );

    return res;
}
