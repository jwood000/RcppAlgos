#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsDesign.h"
#include "Partitions/PartitionsCount.h"
#include "ComboGroup/ComboGroupClass.h"
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

    bool IsRep = CppConvert::convertFlag(RisRep, "repetition");
    bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_class computedRowsMpz;

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    return CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
}

[[cpp11::register]]
SEXP PartitionsCount(SEXP Rtarget, SEXP Rv, SEXP Rm,
                     SEXP RisRep, SEXP RFreqs, SEXP RcompFun,
                     SEXP Rlow, SEXP Rtolerance,
                     SEXP RPartDesign, SEXP Rshow,
                     SEXP RIsComposition, SEXP RIsWeak) {
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
    bool IsRep = CppConvert::convertFlag(RisRep, "repetition");
    const bool bDesign = CppConvert::convertFlag(RPartDesign,
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

    part.isRep   = IsRep;
    part.isMult  = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    part.isWeak  = CppConvert::convertFlag(RIsWeak, "weak");
    part.isComp  = CppConvert::convertFlag(RIsComposition, "composition");
    part.isComb  = !part.isComp;

    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals,
                    funDbl, part, ctype, n, m, compVec, mainFun,
                    mainFun, myType, Rtarget, RcompFun, Rtolerance,
                    Rlow, true);

    if (!part.numUnknown) {
        if (bDesign) {
            bool Verbose = CppConvert::convertFlag(Rshow, "showDetail");
            return GetDesign(part, ctype, n, Verbose);
        } else {
            return CppConvert::GetCount(part.isGmp, part.bigCount,
                                          part.count);
        }
    } else if (bDesign) {
        cpp11::stop("No design available for this case!");
    } else {
        cpp11::stop("The count is unknown for this case.\n To get the"
                    " total number, generate all results!");
    }
}

[[cpp11::register]]
SEXP ComboGroupsCountCpp(SEXP Rv, SEXP RNumGroups, SEXP RGrpSize) {

    int n;
    std::vector<int> vInt;
    std::vector<double> vNum;
    VecType myType = VecType::Integer;

    std::unique_ptr<ComboGroup> CmbGrpCls = GroupPrep(
        vInt, vNum, n, myType, Rv, RNumGroups, RGrpSize
    );

    CmbGrpCls->SetCount();
    return CmbGrpCls->GetCount();
}
