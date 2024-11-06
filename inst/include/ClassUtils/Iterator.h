#pragma once

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

public:

    Iterator(SEXP Rv, VecType typePass, SEXP RcompRow, int RmaxThreads,
             SEXP RnThreads, bool Rparallel, bool IsGmp);

    virtual ~Iterator() = default;
    virtual void startOver() = 0;
    virtual SEXP nextComb() = 0;
    virtual SEXP prevComb() = 0;
    virtual SEXP nextNumCombs(SEXP RNum) = 0;
    virtual SEXP prevNumCombs(SEXP RNum) = 0;
    virtual SEXP nextGather() = 0;
    virtual SEXP prevGather() = 0;
    virtual SEXP currComb() = 0;
    virtual SEXP randomAccess(SEXP RindexVec) = 0;
    virtual SEXP front() = 0;
    virtual SEXP back() = 0;
    virtual SEXP summary() = 0;
    SEXP sourceVector() const;
};
