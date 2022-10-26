#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/RankPartition.h"
#include "RankUtils.h"

void RankPartsResults(std::vector<mpz_class> &bigRes, int* intRes, double* dblRes,
                      std::vector<int> &idx, rankPartsPtr rankFun,
                      int tar, int m, int cap, int strtLen,
                      int numResults, bool IsGmp, bool IsInteger) {

    mpz_class mpzIdx(0);

    if (IsGmp) {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
            ++mpzIdx;
            bigRes[i] = mpzIdx;
        }
    } else {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
            ++dblIdx;
            if (IsInteger) intRes[i] = dblIdx; else dblRes[i] = dblIdx;
        }
    }
}

[[cpp11::register]]
SEXP RankPartitionMain(SEXP RIdx, SEXP Rv, SEXP RisRep,
                       SEXP RFreqs, SEXP Rm, SEXP RcompFun,
                       SEXP Rtarget, SEXP Rtolerance,
                       SEXP RIsComposition, SEXP RIsWeak) {

    int n = 0;
    int m = 0;
    VecType myType = VecType::Integer;

    bool IsRep  = CppConvert::convertFlag(RisRep, "repetition");
    bool IsMult = false;

    std::vector<int> idx;
    std::vector<int> freqs;
    std::vector<int> myReps;
    std::vector<int> vInt;
    std::vector<double> vNum;

    const std::string mainFun = "sum";
    const bool IsComposition  = CppConvert::convertFlag(RIsComposition,
                                                          "IsComposition");

    SetUpRank(RIdx, Rv, RisRep, RFreqs, Rm, idx, freqs, myReps,
              myType, n, m, !IsComposition, IsMult, IsRep);
    SetBasic(Rv, vNum, vInt, n, myType);

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep   = IsRep;
    part.isMult  = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    part.isComp  = IsComposition;
    part.isComb  = !part.isComp;
    part.isWeak  = CppConvert::convertFlag(RIsWeak, "weak");

    cpp11::sexp Rlow = R_NilValue;
    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals, funDbl,
                    part, ctype, n, m, compVec, mainFun, mainFun, myType,
                    Rtarget, RcompFun, Rtolerance, Rlow);

    if (part.ptype == PartitionType::Multiset ||
        part.ptype == PartitionType::CoarseGrained ||
        part.ptype == PartitionType::NotPartition) {

        cpp11::stop("Partition ranking not available for this case.");
    }

    if (part.isComp && !part.isWeak && part.includeZero) {
        for (auto it = idx.rbegin(); it != idx.rend();) {

            bool zero_found = false;

            for (int i = 0; i < m; ++i, ++it) {
                if (vNum[*it] != 0 && zero_found) {
                    cpp11::stop("Malformed composition. If weak = FALSE,"
                                " zero cannot come after nonzero values!");
                }

                if (vNum[*it] == 0) {
                    zero_found = true;
                }
            }
        }
    }

    // See comment in SamplePartitions.cpp
    if (part.numUnknown) PartitionsCount(myReps, part, n, true);
    const int numResults  = Rf_length(RIdx) / m;
    const int bigSampSize = (part.isGmp) ? numResults : 1;
    std::vector<mpz_class> myVec(bigSampSize);

    const int cap     = n - static_cast<int>(part.includeZero);
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const bool IsInt = (part.count <= std::numeric_limits<int>::max());
    cpp11::sexp res_std_int = Rf_allocVector(INTSXP, IsInt ? numResults : 0);
    int* res_int = INTEGER(res_std_int);

    cpp11::sexp res_std_dbl = Rf_allocVector(
        REALSXP, (!IsInt && !part.isGmp) ? numResults : 0
    );
    double* res_dbl = REAL(res_std_dbl);

    if (m == 1) return Rf_ScalarInteger(1);
    const rankPartsPtr rankFun = GetRankPartsFunc(
        part.ptype, part.isGmp, part.isComp
    );

    RankPartsResults(myVec, res_int, res_dbl, idx,
                     rankFun, part.mapTar, m, cap, strtLen,
                     numResults, part.isGmp, IsInt);

    if (IsInt) {
        return res_std_int;
    } else if (part.isGmp) {
        return MpzReturn(myVec, numResults);
    } else {
        return res_std_dbl;
    }
}
