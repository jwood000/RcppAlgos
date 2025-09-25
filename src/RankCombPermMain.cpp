#include "ComputedCount.h"
#include "RankResult.h"
#include "RankUtils.h"
#include <thread>

template <typename T>
void RankCPResultsGeneric(
    std::vector<T> &res, std::vector<int> &idx, const std::vector<int> &myReps,
    rankResultPtr rankFun, int m, int lenV, int numResults
) {

    mpz_class mpzIdx(0);

    for (int i = 0, j = 0; i < numResults; ++i, j += m) {
        double dblIdx = 0;
        rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
        res[i] = RankResultTraits<T>::convert(dblIdx, mpzIdx);
    }
}

template <typename T>
void RankResults(
    T* res, std::vector<int> &idx, const std::vector<int> &myReps,
    rankResultPtr rankFun, int m, int lenV, int numResults
) {

    mpz_class mpzIdx(0);

    for (int i = 0, j = 0; i < numResults; ++i, j += m) {
        double dblIdx = 0;
        rankFun(idx.begin() + j, lenV, m, dblIdx, mpzIdx, myReps);
        ++dblIdx;
        res[i] = dblIdx;
    }
}

template <typename T>
void RankCPParallel(
    std::vector<std::vector<T>> &res_sec, std::vector<T> &last_res,
    std::vector<int> &idx, const std::vector<int> &myReps,
    rankResultPtr rankFun, int m, int lenV, int numResults,
    int nThreads, int stepSize, int last_sec
) {

    std::vector<std::thread> threads;

    std::vector<std::vector<int>> idx_section(nThreads);
    Create2D(idx, idx_section, stepSize, m, nThreads);

    for (int j = 0; j < (nThreads - 1); ++j) {
        threads.emplace_back(
            std::cref(RankCPResultsGeneric<T>), std::ref(res_sec[j]),
            std::ref(idx_section[j]), std::cref(myReps),
            rankFun, m, lenV, stepSize
        );
    }

    threads.emplace_back(
        std::cref(RankCPResultsGeneric<T>), std::ref(last_res),
        std::ref(idx_section.back()), std::cref(myReps),
        rankFun, m, lenV, last_sec
    );

    for (auto& thr: threads) {
        thr.join();
    }
}

template <typename T>
void RankCPStd(
    T* res, std::vector<int> &idx, const std::vector<int> &myReps,
    rankResultPtr rankFun, int m, int lenV, int numResults, int nThreads
) {

    if (nThreads > 1) {
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<T>> res_sec(
            (nThreads - 1), std::vector<T>(stepSize)
        );

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<T> last_res(last_sec);

        RankCPParallel(res_sec, last_res, idx, myReps, rankFun, m,
                       lenV, numResults, nThreads, stepSize, last_sec);

        int offset = 0;

        for (int j = 0; j < (nThreads - 1); ++j) {
            std::copy(res_sec[j].begin(), res_sec[j].end(), res + offset);
            offset += stepSize;
        }

        std::copy(last_res.begin(), last_res.end(), res + offset);
    } else {
        RankResults(res, idx, myReps, rankFun, m, lenV, numResults);
    }
}

void RankCPGmp(
    std::vector<mpz_class> &res, std::vector<int> &idx,
    const std::vector<int> &myReps, rankResultPtr rankFun,
    int m, int lenV, int numResults, int nThreads
) {

    if (nThreads > 1) {
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<mpz_class>> res_sec(
            (nThreads - 1), std::vector<mpz_class>(stepSize)
        );

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<mpz_class> last_res(last_sec);

        RankCPParallel(res_sec, last_res, idx, myReps, rankFun, m,
                       lenV, numResults, nThreads, stepSize, last_sec);

        int offset = 0;

        for (int j = 0; j < (nThreads - 1); ++j) {
            std::copy(res_sec[j].begin(), res_sec[j].end(),
                      res.begin() + offset);
            offset += stepSize;
        }

        std::copy(last_res.begin(), last_res.end(), res.begin() + offset);
    } else {
        RankCPResultsGeneric(res, idx, myReps, rankFun, m, lenV, numResults);
    }
}

[[cpp11::register]]
SEXP RankCombPerm(SEXP RIdx, SEXP Rv, SEXP RisRep,
                  SEXP RFreqs, SEXP Rm, SEXP RIsComb,
                  SEXP RNumThreads, SEXP RmaxThreads) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

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
    const int numResults = Rf_length(RIdx) / m;

    const int limit = 2;
    bool DummyPar = false; // Not used. See notes in SamplePartitions.cpp
    SetThreads(DummyPar, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    const rankResultPtr rankFun = GetRankResultFunc(
        IsComb, IsMult, IsRep, IsGmp
    );

    if (IsGmp) {
        std::vector<mpz_class> myVec(numResults);
        RankCPGmp(myVec, idx, myReps, rankFun, m, n, numResults, nThreads);
        return MpzReturn(myVec, numResults);
    } else if (computedRows <= std::numeric_limits<int>::max()) {
        cpp11::sexp res_std_int = Rf_allocVector(INTSXP, numResults);
        int* res_int = INTEGER(res_std_int);
        RankCPStd(res_int, idx, myReps, rankFun, m, n, numResults, nThreads);
        return res_std_int;
    } else {
        cpp11::sexp res_std_dbl = Rf_allocVector(REALSXP, numResults);
        double* res_dbl = REAL(res_std_dbl);
        RankCPStd(res_dbl, idx, myReps, rankFun, m, n, numResults, nThreads);
        return res_std_dbl;
    }
}
