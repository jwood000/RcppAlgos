#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/RankPartition.h"
#include "RankUtils.h"
#include <thread>

template <typename T>
void RankPartsResultsGeneric(
    std::vector<T>& res, std::vector<int>& idx, rankPartsPtr rankFun,
    int tar, int m, int cap, int strtLen, int numResults
) {

    mpz_class mpzIdx(0);

    for (int i = 0, j = 0; i < numResults; ++i, j += m) {
        double dblIdx = 0;
        rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
        res[i] = RankResultTraits<T>::convert(dblIdx, mpzIdx);
    }
}

template <typename T>
void RankPartsResults(
    T* res, std::vector<int> &idx, rankPartsPtr rankFun,
    int tar, int m, int cap, int strtLen, int numResults
) {

    mpz_class mpzIdx(0);

    for (int i = 0, j = 0; i < numResults; ++i, j += m) {
        double dblIdx = 0;
        rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
        res[i] = (dblIdx + 1);
    }
}

template <typename T>
void RankPartsParallel(
    std::vector<std::vector<T>> &res_sec, std::vector<T> &last_res,
    std::vector<int> &idx, rankPartsPtr rankFun, int tar, int m, int cap,
    int strtLen, int numResults, int nThreads, int stepSize, int last_sec
) {

    std::vector<std::thread> threads;

    std::vector<std::vector<int>> idx_section(nThreads);
    Create2D(idx, idx_section, stepSize, m, nThreads);

    for (int j = 0; j < (nThreads - 1); ++j) {
        threads.emplace_back(
            std::cref(RankPartsResultsGeneric<T>), std::ref(res_sec[j]),
            std::ref(idx_section[j]), rankFun, tar,
            m, cap, strtLen, stepSize
        );
    }

    threads.emplace_back(
        std::cref(RankPartsResultsGeneric<T>), std::ref(last_res),
        std::ref(idx_section.back()), rankFun, tar, m, cap,
        strtLen, last_sec
    );

    for (auto& thr: threads) {
        thr.join();
    }
}

template <typename T>
void RankPartsStd(
    T* res, std::vector<int> &idx, rankPartsPtr rankFun, int tar,
    int m, int cap, int strtLen, int numResults, int nThreads
) {

    if (nThreads > 1) {
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<T>> res_sec(
            (nThreads - 1), std::vector<T>(stepSize)
        );

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<T> last_res(last_sec);

        RankPartsParallel(res_sec, last_res, idx, rankFun, tar, m, cap,
                          strtLen, numResults, nThreads, stepSize, last_sec);

        int offset = 0;

        for (int j = 0; j < (nThreads - 1); ++j) {
            std::copy(res_sec[j].begin(), res_sec[j].end(), res + offset);
            offset += stepSize;
        }

        std::copy(last_res.begin(), last_res.end(), res + offset);
    } else {
        RankPartsResults(res, idx, rankFun, tar,
                         m, cap, strtLen, numResults);
    }
}

void RankPartsGmp(
    std::vector<mpz_class> &res, std::vector<int> &idx,
    rankPartsPtr rankFun, int tar, int m, int cap,
    int strtLen, int numResults, int nThreads
) {

    if (nThreads > 1) {
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<mpz_class>> res_sec(
            (nThreads - 1), std::vector<mpz_class>(stepSize)
        );

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<mpz_class> last_res(last_sec);

        RankPartsParallel(res_sec, last_res, idx, rankFun, tar, m, cap,
                          strtLen, numResults, nThreads, stepSize, last_sec);

        int offset = 0;

        for (int j = 0; j < (nThreads - 1); ++j) {
            std::copy(res_sec[j].begin(), res_sec[j].end(),
                      res.begin() + offset);
            offset += stepSize;
        }

        std::copy(last_res.begin(), last_res.end(), res.begin() + offset);
    } else {
        RankPartsResultsGeneric(res, idx, rankFun, tar,
                                m, cap, strtLen, numResults);
    }
}

[[cpp11::register]]
SEXP RankPartitionMain(SEXP RIdx, SEXP Rv, SEXP RisRep,
                       SEXP RFreqs, SEXP Rm, SEXP RcompFun,
                       SEXP Rtarget, SEXP Rtolerance,
                       SEXP RIsComposition, SEXP RIsWeak,
                       SEXP RNumThreads, SEXP RmaxThreads) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    bool IsRep  = CppConvert::convertFlag(RisRep, "repetition");
    bool IsMult = false;

    std::vector<int> idx;
    std::vector<int> freqs;
    std::vector<int> myReps;
    std::vector<int> vInt;
    std::vector<double> vNum;

    const std::string mainFun = "sum";
    bool IsComp = CppConvert::convertFlag(RIsComposition, "IsComposition");

    SetUpRank(RIdx, Rv, RisRep, RFreqs, Rm, idx, freqs, myReps,
              myType, n, m, !IsComp, IsMult, IsRep);
    SetBasic(Rv, vNum, vInt, n, myType);

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;
    InitialSetupPartDesign(part, RIsWeak, RIsComposition, IsRep,
                           IsMult, Rf_isNull(Rm), !IsComp);

    cpp11::sexp Rlow = R_NilValue;
    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals, funDbl,
                    part, ctype, n, m, compVec, mainFun, mainFun, myType,
                    Rtarget, RcompFun, Rtolerance, Rlow);

    if (part.ptype == PartitionType::CoarseGrained ||
        part.ptype == PartitionType::NotPartition  ||
        part.ptype == PartitionType::NoSolution    ||
        part.ptype == PartitionType::Multiset) {

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

    const int cap     = n - static_cast<int>(part.includeZero);
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const int limit = 2;
    bool DummyPar = false; // Not used. See notes in SamplePartitions.cpp
    SetThreads(DummyPar, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    const rankPartsPtr rankFun = GetRankPartsFunc(part.ptype, part.isGmp);

    if (part.isGmp) {
        std::vector<mpz_class> myVec(numResults);
        RankPartsGmp(myVec, idx, rankFun, part.mapTar, m,
                     cap, strtLen, numResults, nThreads);
        return MpzReturn(myVec, numResults);
    } else if (part.count <= std::numeric_limits<int>::max()) {
        cpp11::sexp res_std_int = Rf_allocVector(INTSXP, numResults);
        int* res_int = INTEGER(res_std_int);

        RankPartsStd(res_int, idx, rankFun, part.mapTar,
                     m, cap, strtLen, numResults, nThreads);

        return res_std_int;
    } else {
        cpp11::sexp res_std_dbl = Rf_allocVector(REALSXP, numResults);
        double* res_dbl = REAL(res_std_dbl);

        RankPartsStd(res_dbl, idx, rankFun, part.mapTar,
                     m, cap, strtLen, numResults, nThreads);

        return res_std_dbl;
    }
}
