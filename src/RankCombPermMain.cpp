#include "ComputedCount.h"
#include "RankResult.h"
#include "RankUtils.h"

void RankResults(std::vector<mpz_class> &bigRes, int* intRes, double* dblRes,
                 std::vector<int> &idx, const std::vector<int> &myReps,
                 rankResultPtr rankFun, int m, int lenV, int numResults,
                 bool IsGmp, bool IsInteger) {

    mpz_class mpzIdx(0);

    if (IsGmp) {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
            ++mpzIdx;
            bigRes[i] = mpzIdx;
        }
    } else {
        for (int i = 0, j = 0; i < numResults; ++i, j += m) {
            double dblIdx = 0;
            rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
            ++dblIdx;
            if (IsInteger) intRes[i] = dblIdx; else dblRes[i] = dblIdx;
        }
    }
}

[[cpp11::register]]
SEXP RankCombPerm(SEXP RIdx, SEXP Rv, SEXP RisRep,
                  SEXP RFreqs, SEXP Rm, SEXP RIsComb) {

    int n = 0;
    int m = 0;
    VecType myType = VecType::Integer;

    bool IsRep = CppConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");
    bool IsMult = false;

    std::vector<int> idx;
    std::vector<int> freqs;
    std::vector<int> myReps;

    SetUpRank(RIdx, Rv, RisRep, RFreqs, Rm, idx, freqs,
              myReps, myType, n, m, IsComb, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp     = (computedRows > Significand53);
    const bool IsInteger = (computedRows <= std::numeric_limits<int>::max());

    const int numResults = Rf_length(RIdx) / m;
    const rankResultPtr rankFun = GetRankResultFunc(
        IsComb, IsMult, IsRep, IsGmp
    );
    cpp11::sexp res_std_int = Rf_allocVector(INTSXP, (IsInteger) ?
                                                 numResults : 0);
    int* res_int = INTEGER(res_std_int);

    cpp11::sexp res_std_dbl = Rf_allocVector(
        REALSXP, (!IsInteger && !IsGmp) ? numResults : 0
    );
    double* res_dbl = REAL(res_std_dbl);

    const int bigNumResults = IsGmp ? numResults : 0;
    std::vector<mpz_class> myVec(bigNumResults);

    RankResults(myVec, res_int, res_dbl, idx, myReps,
                rankFun, m, n, numResults, IsGmp, IsInteger);

    if (IsInteger) {
        return res_std_int;
    } else if (IsGmp) {
        return MpzReturn(myVec, numResults);
    } else {
        return res_std_dbl;
    }
}
