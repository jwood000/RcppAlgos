#pragma once

#include "SetUpUtils.h"

template <typename T>
struct RankResultTraits;

template <>
struct RankResultTraits<int> {
    static int convert(double dblIdx, const mpz_class&) {
        return static_cast<int>(dblIdx + 1);
    }
};

template <>
struct RankResultTraits<double> {
    static double convert(double dblIdx, const mpz_class&) {
        return dblIdx + 1;
    }
};

template <>
struct RankResultTraits<mpz_class> {
    static mpz_class convert(double, const mpz_class& mpzIdx) {
        return mpzIdx + 1;
    }
};

template <typename T>
void Create2D(const std::vector<int> &idx,
              std::vector<std::vector<T>> &v, int stepSize,
              int m, int nThreads);

void SetUpRank(SEXP RIdx, SEXP Rv, SEXP RisRep, SEXP RFreqs, SEXP Rm,
               std::vector<int> &idx, std::vector<int> &freqs,
               std::vector<int> &myReps, VecType &myType, int &n, int &m,
               bool IsComb, bool &IsMult, bool &IsRep);

SEXP MpzReturn(const std::vector<mpz_class> &myVec, int numResults);
