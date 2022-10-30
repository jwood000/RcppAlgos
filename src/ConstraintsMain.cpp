#include "Constraints/GetContraints.h"
#include "Partitions/NthPartition.h"
#include "ComputedCount.h"
#include "CheckReturn.h"

[[cpp11::register]]
SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                         SEXP Rlow, SEXP Rhigh, SEXP RmainFun,
                         SEXP RcompFun, SEXP Rtarget, SEXP RIsComb,
                         SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads,
                         SEXP RmaxThreads, SEXP Rtolerance,
                         SEXP RIsComposition, SEXP RIsWeak) {
    int n = 0;
    int m = 0;
    int nRows = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool KeepRes  = CppConvert::convertFlag(RKeepRes, "keepResults");
    bool Parallel = CppConvert::convertFlag(Rparallel, "Parallel");
    bool IsRep    = CppConvert::convertFlag(RisRep, "repetition");

    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");
    const bool IsConstrained = CheckConstrnd(RmainFun, RcompFun, Rtarget);

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, IsConstrained);

    if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
        cpp11::stop("contraintFun must be one of the following:"
                 " 'prod', 'sum', 'mean', 'max', or 'min'");
    }

    const std::string funTest(CHAR(STRING_ELT(RmainFun, 0)));
    const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), funTest);

    if (funIt == mainFunSet.end()) {
        cpp11::stop("contraintFun must be one of the following:"
                 " 'prod', 'sum', 'mean', 'max', or 'min'");
    }

    // Must be defined inside IsInteger check as tarVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> tarIntVals;
    const std::string mainFun = funTest == "mean" ? "sum" : funTest;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> tarVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep   = IsRep;
    part.isMult  = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    part.isWeak  = CppConvert::convertFlag(RIsWeak, "weak");
    part.isComp  = CppConvert::convertFlag(RIsComposition,
                                             "IsComposition");
    part.isComb = IsComb;

    if (IsConstrained) {
        ConstraintSetup(vNum, myReps, tarVals, vInt, tarIntVals,
                        funDbl, part, ctype, n, m, compVec, mainFun,
                        funTest, myType, Rtarget, RcompFun,
                        Rtolerance, Rlow);
    }

    const bool usePartCount = part.isPart &&
                              !part.isGmp &&
                              !part.numUnknown;

    const double computedRows = usePartCount ? part.count :
        GetComputedRows(IsMult, IsComb, IsRep, n, m, Rm, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_class computedRowsMpz;

    if (IsGmp && part.isPart) {
        computedRowsMpz = part.bigCount;
    } else if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    // This variable is used in determining the number of results. If the
    // output is constrained and the ConstraintType is "General" or
    // "PartitionEsque", it means we really don't know how many results
    // we have. The computedRows above is a strict upper bound but not
    // necessarily the least upper bound. In these cases, we don't want
    // to unnecessarily throw an error when computedRows exceeds 2^31 - 1.
    const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
                            ctype == ConstraintType::SpecialCnstrnt ||
                            ctype == ConstraintType::General        ||
                            (part.isPart && part.numUnknown);

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz, upperMpz, computedRowsMpz, computedRows);

    std::vector<int> startZ(m);
    const int cap     = n - part.includeZero;
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    if (ctype < ConstraintType::PartMapping) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz, IsRep, IsMult, IsGmp);
    } else {
        if (bLower) {
            const nthPartsPtr nthPartFun = GetNthPartsFunc(
                part.ptype, IsGmp, part.isComp
            );
            startZ = nthPartFun(part.mapTar, part.width, cap,
                                strtLen, lower, lowerMpz);

            if (ctype == ConstraintType::PartStandard && !part.includeZero) {
                for (auto &z_i: startZ) {
                    ++z_i;
                }
            }
        } else {
            startZ = part.startZ;
        }
    }

    // This is used when we are unable to calculate the number of results
    // upfront (E.g. comboGeneral(rnorm(10), 5, constraintFun = "sum,
    //                            comparisonFun = "<=", limitConstraints = 1))
    double userNum = 0;
    const bool bSetNum = !numUnknown ||
        ctype == ConstraintType::SpecialCnstrnt;

    SetNumResults(IsGmp, bLower, bUpper, bSetNum, upperMpz,
                  lowerMpz, lower, upper, computedRows,
                  computedRowsMpz, nRows, userNum);

    int nThreads   = 1;
    int maxThreads = 1;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    const int limit = (part.isPart) ?
    ((part.ptype == PartitionType::RepCapped   ||
      part.ptype == PartitionType::DstctCapped ||
      part.ptype == PartitionType::DstctCappedMZ) ? 150000 : 40000) : 20000;

    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RnThreads, limit);

    cpp11::sexp res = GetConstraints(
      part, compVec, freqs, myReps, vNum, vInt, tarVals, tarIntVals,
      startZ, mainFun, funTest, funDbl, lower, lowerMpz, userNum,
      ctype, myType, nThreads, nRows, n, strtLen, cap, m, IsComb,
      Parallel, IsGmp, IsRep, IsMult, bUpper, KeepRes, numUnknown
    );

    return res;
}
