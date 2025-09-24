#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/RankPartition.h"
#include "RankUtils.h"
#include <thread>

template <typename T>
struct RankResultTraits;

// int and double: same formula
template <>
struct RankResultTraits<int> {
    using result_type = int;
    static result_type convert(double dblIdx, const mpz_class&) {
        return static_cast<int>(dblIdx + 1);
    }
};

template <>
struct RankResultTraits<double> {
    using result_type = double;
    static result_type convert(double dblIdx, const mpz_class&) {
        return dblIdx + 1;
    }
};

// mpz_class: use the mpz_class result
template <>
struct RankResultTraits<mpz_class> {
    using result_type = mpz_class;
    static result_type convert(double, const mpz_class& mpzIdx) {
        return mpzIdx + 1;
    }
};

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

void RankPartsResultsGmp(
    std::vector<mpz_class> &res, std::vector<int> &idx, rankPartsPtr rankFun,
    int tar, int m, int cap, int strtLen, int numResults
) {

    mpz_class mpzIdx(0);

    for (int i = 0, j = 0; i < numResults; ++i, j += m) {
        double dblIdx = 0;
        rankFun(idx.begin() + j, tar, m, cap, strtLen, dblIdx, mpzIdx);
        res[i] = (mpzIdx + 1);
    }
}

template <typename T>
void RankPartsWorker(
    std::vector<T> &res, std::vector<int> &idx, rankPartsPtr rankFun,
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
void Create2D(const std::vector<int> &idx,
              std::vector<std::vector<T>> &v, int stepSize,
              int m, int nThreads) {

    int strt = 0;
    const int sectionWidth = stepSize * m;

    for (int i = 0; i < (nThreads - 1); ++i) {
        v[i].assign(idx.begin() + strt, idx.begin() + strt + sectionWidth);
        strt += sectionWidth;
    }

    v.back().assign(idx.begin() + strt, idx.end());
}

template <typename T>
void RankPartsStdMain(
    T* res, std::vector<int> &idx, rankPartsPtr rankFun, int tar,
    int m, int cap, int strtLen, int numResults, int nThreads
) {

    if (nThreads > 1) {
        std::vector<std::thread> threads;
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<T>> res_sec(
            (nThreads - 1), std::vector<T>(stepSize)
        );

        std::vector<std::vector<int>> idx_section(nThreads);
        Create2D(idx, idx_section, stepSize, m, nThreads);

        for (int j = 0; j < (nThreads - 1); ++j) {
            threads.emplace_back(
                std::cref(RankPartsWorker<T>), std::ref(res_sec[j]),
                std::ref(idx_section[j]), rankFun, tar,
                m, cap, strtLen, stepSize
            );
        }

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<T> last_res(last_sec);

        threads.emplace_back(
            std::cref(RankPartsWorker<T>), std::ref(last_res),
            std::ref(idx_section.back()), rankFun, tar, m, cap,
            strtLen, last_sec
        );

        for (auto& thr: threads) {
            thr.join();
        }

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

SEXP RankPartsStdGmp(
    std::vector<int> &idx, rankPartsPtr rankFun, int tar,
    int m, int cap, int strtLen, int numResults, int nThreads
) {

    if (nThreads > 1) {
        std::vector<std::thread> threads;
        const int stepSize = numResults / nThreads;

        std::vector<std::vector<mpz_class>> res_sec(
            (nThreads - 1), std::vector<mpz_class>(stepSize)
        );

        std::vector<std::vector<int>> idx_section(nThreads);
        Create2D(idx, idx_section, stepSize, m, nThreads);

        for (int j = 0; j < (nThreads - 1); ++j) {
            threads.emplace_back(
                std::cref(RankPartsResultsGmp), std::ref(res_sec[j]),
                std::ref(idx_section[j]), rankFun, tar,
                m, cap, strtLen, stepSize
            );
        }

        const int last_sec = numResults - ((nThreads - 1) * stepSize);
        std::vector<mpz_class> last_res(last_sec);

        threads.emplace_back(
            std::cref(RankPartsResultsGmp), std::ref(last_res),
            std::ref(idx_section.back()), rankFun, tar, m, cap,
            strtLen, last_sec
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        std::vector<mpz_class> myVec(numResults);
        RankPartsResultsGmp(myVec, idx, rankFun, tar,
                            m, cap, strtLen, numResults);
        return MpzReturn(myVec, numResults);
    }

    return Rf_ScalarInteger(1);
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

    bool Parallel = false; // This will be set in SetThreads below. For the
    // partition and composition functions we don't have a Parallel
    // argument. The goal is to eventually phase out this argument.

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

    const int cap     = n - static_cast<int>(part.includeZero);
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const int limit = 2;
    SetThreads(Parallel, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    const rankPartsPtr rankFun = GetRankPartsFunc(part.ptype, part.isGmp);

    if (part.isGmp) {
        return RankPartsStdGmp(idx, rankFun, part.mapTar, m,
                               cap, strtLen, numResults, nThreads);
    } else if (part.count <= std::numeric_limits<int>::max()) {
        cpp11::sexp res_std_int = Rf_allocVector(INTSXP, numResults);
        int* res_int = INTEGER(res_std_int);

        RankPartsStdMain(
            res_int, idx, rankFun, part.mapTar,
            m, cap, strtLen, numResults, nThreads
        );

        return res_std_int;
    } else {
        cpp11::sexp res_std_dbl = Rf_allocVector(REALSXP, numResults);
        double* res_dbl = REAL(res_std_dbl);

        RankPartsStdMain(
            res_dbl, idx, rankFun, part.mapTar,
            m, cap, strtLen, numResults, nThreads
        );

        return res_std_dbl;
    }
}
