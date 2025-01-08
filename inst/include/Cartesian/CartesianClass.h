#pragma once

#include "Cartesian/GetProduct.h"
#include "ClassUtils/Iterator.h"

class CartesianClass : public Iterator {
private:

    const cpp11::list RList;

    const std::vector<int> idx;
    const std::vector<int> lenGrps;
    const std::vector<int> typeCheck;
    const std::vector<int> IsFactor;

    const std::vector<int> intVec;
    const std::vector<double> dblVec;
    const std::vector<int> boolVec;
    const std::vector<Rcomplex> cmplxVec;
    const std::vector<Rbyte> rawVec;
    const cpp11::strings charVec;

    const bool IsDF;
    const int nCols;

    std::vector<int> z;
    std::vector<int> lenNxtPr;
    const VecType myType;

    SEXP VectorReturn();
    SEXP SingleReturn();
    SEXP GeneralReturn(int numResults);

public:

    CartesianClass(
        SEXP Rv_RList, SEXP RcompRows, int RmaxThreads, SEXP RnumThreads,
        bool Rparallel, bool RIsGmp, const std::vector<int> &Ridx,
        const std::vector<int> &RtypeCheck, const std::vector<int> &RIsFactor,
        const std::vector<int> &RintVec, const std::vector<double> &RdblVec,
        const std::vector<int> &RlglVec, const std::vector<Rcomplex> &RcplxVec,
        const std::vector<Rbyte> &RrawVec, const cpp11::strings &RcharVec,
        const std::vector<int> &RlenGrps,
        bool RisDF, int RnCols, VecType RmyType
    );

    void startOver();
    SEXP nextIter();
    SEXP nextNumIters(SEXP RNum);
    SEXP nextGather();
    SEXP currIter();
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP summary();

    // Not currently implemented for this class.
    SEXP prevIter();
    SEXP prevNumIters(SEXP RNum);
    SEXP prevGather();
};
