#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/RankPartition.h"
#include "Cpp14MakeUnique.h"
#include "RankUtils.h"

void RankPartsResults(mpz_t* bigRes, int* intRes, double* dblRes,
                      std::vector<int> &idx, rankPartsPtr rankFun,
                      int tar, int m, int cap, int strtLen,
                      int numResults, bool IsGmp, bool IsInteger) {

    mpz_t mpzIdx;
    mpz_init_set_ui(mpzIdx, 0u);

    if (IsGmp) {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
            mpz_add_ui(mpzIdx, mpzIdx, 1u);
            mpz_set(bigRes[i], mpzIdx);
        }
    } else {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
            ++dblIdx;
            if (IsInteger) intRes[i] = dblIdx; else dblRes[i] = dblIdx;
        }
    }

    mpz_clear(mpzIdx);
}

[[cpp11::register]]
SEXP RankPartitionMain(SEXP RIdx, SEXP Rv, SEXP RisRep,
                       SEXP RFreqs, SEXP Rm, SEXP RcompFun,
                       SEXP Rtarget, SEXP Rtolerance) {

    int n = 0;
    int m = 0;
    VecType myType = VecType::Integer;

    bool IsRep  = CleanConvert::convertFlag(RisRep, "repetition");
    bool IsMult = false;

    std::vector<int> idx;
    std::vector<int> freqs;
    std::vector<int> myReps;
    std::vector<int> vInt;
    std::vector<double> vNum;

    SetUpRank(RIdx, Rv, RisRep, RFreqs, Rm, idx, freqs,
              myReps, myType, n, m, true, IsMult, IsRep);
    SetBasic(Rv, vNum, vInt, n, myType);

    const std::string mainFun = "sum";

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep       = IsRep;
    part.isMult      = IsMult;
    part.mIsNull     = Rf_isNull(Rm);
    cpp11::sexp Rlow = R_NilValue;

    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals,
                    funDbl, part, ctype, n, m, compVec, mainFun, mainFun,
                    myType, Rtarget, RcompFun, Rtolerance, Rlow, true, false);

    if (part.ptype == PartitionType::Multiset ||
        part.ptype == PartitionType::CoarseGrained ||
        part.ptype == PartitionType::NotPartition) {

        cpp11::stop("Partition ranking not available for this case.");
    }

    const int numResults  = Rf_length(RIdx) / m;
    const int bigSampSize = (part.isGmp) ? numResults : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);

    for (int i = 0; i < bigSampSize; ++i) {
        mpz_init(myVec[i]);
    }

    const int cap     = n - static_cast<int>(part.includeZero);
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const bool IsInteger = (part.count <= std::numeric_limits<int>::max());
    cpp11::sexp res_std_int = Rf_allocVector(INTSXP, (IsInteger) ?
                                                 numResults : 0);
    int* res_int = INTEGER(res_std_int);

    cpp11::sexp res_std_dbl = Rf_allocVector(REALSXP, (!IsInteger && !part.isGmp)
                                                 ? numResults : 0);
    double* res_dbl = REAL(res_std_dbl);

    if (m == 1) return Rf_ScalarInteger(1);
    const rankPartsPtr rankFun = GetRankPartsFunc(part.ptype, part.isGmp);

    RankPartsResults(myVec.get(), res_int, res_dbl, idx,
                     rankFun, part.mapTar, m, cap, strtLen,
                     numResults, part.isGmp, IsInteger);

    if (IsInteger) {
        return res_std_int;
    } else if (part.isGmp) {
        return MpzReturn(myVec.get(), numResults);
    } else {
        return res_std_dbl;
    }
}
