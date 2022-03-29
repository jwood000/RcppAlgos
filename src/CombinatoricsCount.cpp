#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsDesign.h"
#include "Partitions/PartitionsCount.h"
#include "ComboGroupsUtils.h"
#include "ImportExportMPZ.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

[[cpp11::register]]
SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep,
                        SEXP RFreqs, SEXP RIsComb) {

    int n = 0;
    int m = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;
    std::vector<double> vNum;

    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    return CleanConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
}

[[cpp11::register]]
SEXP PartitionsCount(SEXP Rtarget, SEXP Rv, SEXP Rm,
                     SEXP RisRep, SEXP RFreqs, SEXP RcompFun,
                     SEXP Rlow, SEXP Rtolerance,
                     SEXP RPartDesign, SEXP Rshow) {
    int n = 0;
    int m = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    const bool IsConstrained = true;
    const std::string mainFun = "sum";
    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    const bool bDesign = CleanConvert::convertFlag(RPartDesign,
                                                   "PartitionsDesign");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, IsConstrained);

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> targetVals;

    ConstraintType ctype;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);

    if (IsConstrained) {
        ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals,
                        funDbl, part, ctype, n, m, compVec, mainFun,
                        mainFun, myType, Rtarget, RcompFun, Rtolerance,
                        Rlow, true, true);
    }

    if (part.ptype != PartitionType::CoarseGrained &&
        part.ptype != PartitionType::NotPartition) {

        if (bDesign) {
            bool Verbose = CleanConvert::convertFlag(Rshow, "showDetail");
            return GetDesign(part, ctype, n, Verbose);
        } else {
            return CleanConvert::GetCount(part.isGmp, part.bigCount,
                                          part.count);
        }
    } else if (bDesign) {
        cpp11::stop("No design available for this case!");
    } else {
        cpp11::stop("The count is unknown for this case. To get the total"
                 " number, generate all results!");
    }
}

[[cpp11::register]]
SEXP ComboGroupsCountCpp(SEXP Rv, SEXP RNumGroups) {

    int n, numGroups;
    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RNumGroups, numGroups,
                                   VecType::Integer, "numGroups");

    std::vector<int> vInt;
    std::vector<double> vNum;

    SetType(myType, Rv);
    SetBasic(Rv, vNum, vInt, n, myType);

    if (n % numGroups != 0) {
        cpp11::stop("The length of v (if v is a vector) or v (if v"
                 " is a scalar) must be divisible by numGroups");
    }

    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > Significand53);

    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);

    if (IsGmp) {
        mpz_set_ui(computedRowMpz, 1);
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    }

    return CleanConvert::GetCount(IsGmp, computedRowMpz, computedRows);
}
