#include "cpp11/list.hpp"

#include "Partitions/BrutePartitionsCountMultiset.h"
#include "ComboGroups/ComboGroupsTemplate.h"
#include "Partitions/MultisetCountClass.h"
#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsDesign.h"
#include "Partitions/PartitionsCount.h"
#include "Cartesian/CartesianUtils.h"
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
    const bool IsGmp = IsBeyondBound(computedRows);
    mpz_class computedRowsMpz;

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    return CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
}

[[cpp11::register]]
SEXP PartitionsCount(
    SEXP Rtarget, SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP RIsComb,
    SEXP RcompFun, SEXP Rlow, SEXP Rtolerance, SEXP RPartDesign, SEXP Rshow,
    SEXP RIsComposition, SEXP RIsWeak
) {

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
    bool bDesign = CppConvert::convertFlag(RPartDesign, "PartitionsDesign");
    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");

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
    InitialSetupPartDesign(part, RIsWeak, RIsComposition, IsRep,
                           IsMult, Rf_isNull(Rm), IsComb);

    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals,
                    funDbl, part, ctype, n, m, compVec, mainFun,
                    mainFun, myType, Rtarget, RcompFun, Rtolerance, Rlow);

    if (!part.numUnknown) {
        if (bDesign) {
            bool Verbose = CppConvert::convertFlag(Rshow, "showDetail");
            return GetDesign(part, ctype, n, Verbose);
        } else {
            return CppConvert::GetCount(part.isGmp, part.bigCount, part.count);
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

    SetType(myType, Rv);
    SetBasic(Rv, R_NilValue, vNum, vInt, n, myType);

    std::unique_ptr<ComboGroupsTemplate> CmbGrpCls =
        GroupPrep(Rv, RNumGroups, RGrpSize, n);

    CmbGrpCls->SetCount();
    return CmbGrpCls->GetCount();
}

[[cpp11::register]]
SEXP ExpandGridCountCpp(cpp11::list RList) {

    const int nCols = Rf_length(RList);
    std::vector<int> lenGrps(nCols);

    for (int i = 0; i < nCols; ++i) {
        lenGrps[i] = Rf_length(RList[i]);
    }

    const double computedRows = CartesianCount(lenGrps);
    const bool IsGmp = IsBeyondBound(computedRows);
    mpz_class computedRowsMpz;

    if (IsGmp) {
        CartesianCountGmp(computedRowsMpz, lenGrps);
    }

    return CppConvert::GetCount(IsGmp, computedRowsMpz, computedRows);
}

[[cpp11::register]]
SEXP PartitionsMultisetCount(
    SEXP Rtarget, SEXP Rv, SEXP Rm, SEXP RFreqs,
    SEXP RIsComb, SEXP RcompFun, SEXP Rlow, SEXP Rtolerance,
    SEXP RIsComposition, SEXP RIsWeak, SEXP RCheckGmp
) {

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
    bool IsRep = false;
    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");
    const bool CheckGmp = CppConvert::convertFlag(RCheckGmp, "checkGmp");

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
    InitialSetupPartDesign(part, RIsWeak, RIsComposition, IsRep,
                           IsMult, Rf_isNull(Rm), IsComb);

    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals,
                    funDbl, part, ctype, n, m, compVec, mainFun,
                    mainFun, myType, Rtarget, RcompFun, Rtolerance, Rlow);

    if (part.isGmp || IsBeyondBound(part.count)) {
        cpp11::stop(
            "This checker only validates counts representable as double."
        );
    }

    if (part.numUnknown) {
        cpp11::stop(
            "The count is unknown for this case.\n"
            "To get the total number, generate all results!"
        );
    }

    if (CheckGmp) {
        const bool vHasZero = (vNum.front() == 0);
        std::vector<int> countReps = myReps;

        if (part.isComp && vHasZero && !part.isWeak) {
            countReps.erase(countReps.begin());
        }

        std::unique_ptr<CountClass> Counter =
            MakeMultisetCount(part.ptype, countReps);

        if (!Counter) {
            cpp11::stop("GMP checker does not support this partition type.");
        }

        if (Counter) {
            std::vector<int> allowed(part.cap);
            std::iota(allowed.begin(), allowed.end(), 1);

            const int strtLen = std::count_if(
                part.startZ.cbegin(), part.startZ.cend(), [](int i){return i > 0;}
            );

            Counter->SetArrSize(part.ptype, part.mapTar, part.width);
            Counter->InitializeMpz();
            Counter->GetCount(
                part.bigCount, part.mapTar, part.width, allowed, strtLen
            );
        }

        cpp11::writable::list res_gmp_lst(3);

        res_gmp_lst[0] = Rf_ScalarReal(part.bigCount.get_d());
        res_gmp_lst[1] = Rf_ScalarReal(part.count);
        res_gmp_lst[2] = Rf_ScalarLogical(part.bigCount.get_d() == part.count);

        res_gmp_lst.names() = {"gmp_algo", "dbl_algo", "check"};
        return res_gmp_lst;
    }

    double result = 0.0;

    switch (part.ptype) {
        case PartitionType::LengthOne:
            result = static_cast<double>(part.solnExist);
            break;

        case PartitionType::Multiset:
            result = part.solnExist ?
                CountPartsMultisetBrute(myReps, part.startZ) : 0.0;
            break;

        case PartitionType::CompMultiset:
            // z[0] is an index referring to the smallest positive value.
            // It must participate in permutation counting.
            result = part.solnExist ?
                CountPartsMultisetBrute(
                    myReps, part.startZ, true, true
                ) : 0.0;
            break;

        case PartitionType::CompMultisetZero:
            // Here index 0 represents an actual zero in v and is padding
            // for non-weak output, so exclude it from permutations.
            result = part.solnExist ?
                CountPartsMultisetBrute(
                    myReps, part.startZ, true, false
                ) : 0.0;
            break;

        case PartitionType::CompMultisetWeak:
        case PartitionType::PrmMultiset:
            // Actual zero, when present, participates in the result.
            result = part.solnExist ?
                CountPartsMultisetBrute(
                    myReps, part.startZ, true, true
                ) : 0.0;
        break;

        default:
            cpp11::stop(
                "Unexpected PartitionType in PartitionsMultisetCount: %s",
                GetPTypeName(part.ptype).c_str()
            );
    }

    cpp11::writable::list res_lst(3);

    res_lst[0] = Rf_ScalarReal(result);
    res_lst[1] = Rf_ScalarReal(part.count);
    res_lst[2] = Rf_ScalarLogical(result == part.count);

    // Assign names to the elements
    res_lst.names() = {"brute", "algo", "check"};
    return res_lst;
}
