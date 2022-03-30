#include "Constraints/GetContraints.h"
#include "Partitions/NthPartition.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "CheckReturn.h"

[[cpp11::register]]
SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                         SEXP Rlow, SEXP Rhigh, SEXP RmainFun,
                         SEXP RcompFun, SEXP Rtarget, SEXP RIsComb,
                         SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads,
                         SEXP RmaxThreads, SEXP Rtolerance) {
    int n = 0;
    int m = 0;
    int nRows = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool KeepRes  = CleanConvert::convertFlag(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertFlag(Rparallel, "Parallel");
    bool IsRep    = CleanConvert::convertFlag(RisRep, "repetition");

    const bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");
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

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);

    if (IsConstrained) {
        ConstraintSetup(vNum, myReps, tarVals, vInt, tarIntVals,
                        funDbl, part, ctype, n, m, compVec, mainFun,
                        funTest, myType, Rtarget, RcompFun,
                        Rtolerance, Rlow, IsComb, false);
    }

    const double computedRows = (part.count > 0) ? part.count :
        GetComputedRows(IsMult, IsComb, IsRep, n, m, Rm, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp && part.isPart) {
        mpz_set(computedRowsMpz, part.bigCount);
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
                            part.numUnknown;

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz,
              computedRows);

    std::vector<int> startZ(m);
    const int cap     = n - part.includeZero;
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    if (ctype < ConstraintType::PartMapping) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz[0], IsRep, IsMult, IsGmp);
    } else {
        if (bLower) {
            nthPartsPtr nthPartFun = GetNthPartsFunc(part.ptype, IsGmp);
            startZ = nthPartFun(part.mapTar, part.width,
                                cap, strtLen, lower, lowerMpz[0]);

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

    SetNumResults(IsGmp, bLower, bUpper, bSetNum, upperMpz[0],
                  lowerMpz[0], lower, upper, computedRows,
                  computedRowsMpz, nRows, userNum);
    mpz_clear(computedRowsMpz);

    int nThreads = 1;
    int maxThreads = 1;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    const int limit = (part.isPart) ?
    ((part.ptype == PartitionType::RepCapped   ||
      part.ptype == PartitionType::DstctCapped ||
      part.ptype == PartitionType::DstctCappedMZ) ? 150000 : 40000) : 20000;

    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RnThreads, limit);

    SEXP res = PROTECT(
        GetConstraints(
          part, compVec, freqs, myReps, vNum, vInt, tarVals, tarIntVals,
          startZ, mainFun, funTest, funDbl, lower, lowerMpz[0], userNum,
          ctype, myType, nThreads, nRows, n, strtLen, cap, m, IsComb,
          Parallel, IsGmp, IsRep, IsMult, bUpper, KeepRes, numUnknown
        )
    );

    mpz_clear(lowerMpz[0]);
    mpz_clear(upperMpz[0]);
    UNPROTECT(1);
    return res;
}
