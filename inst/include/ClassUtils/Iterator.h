#pragma once

#include "SetUpUtils.h"

class Iterator {
protected:

    const SEXP sexpVec;
    const int n;
    int RTYPE;

    const int maxThreads;
    const SEXP sexpNThreads;
    const bool Parallel;

    // IsGmp may change depending on whether we have partitions
    bool IsGmp;
    double computedRows;
    mpz_class computedRowsMpz;

public:

    Iterator(SEXP Rv, SEXP RcompRow, int RmaxThreads,
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
