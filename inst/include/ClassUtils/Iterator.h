#pragma once

#include "ClassUtils/ClassUtils.h"
#include "SetUpUtils.h"

class Iterator {
protected:

    const int n;
    const SEXP sexpVec;
    int RTYPE;
    const VecType myType;

    const int maxThreads;
    const SEXP sexpNThreads;

    const bool Parallel;

    // IsGmp may change depending on whether we have partitions
    bool IsGmp;

    double computedRows;
    mpz_class computedRowsMpz;

    std::vector<int> z;

    double dblTemp;
    mpz_class mpzTemp;

    double dblIndex;
    mpz_class mpzIndex;

    bool prevIterAvailable;
    SEXP ToSeeLast(bool AdjustIdx = true);
    SEXP ToSeeFirst(bool AdjustIdx = true);

public:

    Iterator(SEXP Rv, VecType typePass, SEXP RcompRow, int RmaxThreads,
             SEXP RnThreads, bool Rparallel, bool IsGmp);

    virtual ~Iterator() = default;
    virtual void startOver() = 0;
    virtual SEXP nextIter() = 0;
    virtual SEXP prevIter() = 0;
    virtual SEXP nextNumIters(SEXP RNum) = 0;
    virtual SEXP prevNumIters(SEXP RNum) = 0;
    virtual SEXP nextGather() = 0;
    virtual SEXP prevGather() = 0;
    virtual SEXP currIter() = 0;
    virtual SEXP randomAccess(SEXP RindexVec) = 0;
    virtual SEXP front() = 0;
    virtual SEXP back() = 0;
    virtual SEXP summary() = 0;
    SEXP sourceVector() const;
};
