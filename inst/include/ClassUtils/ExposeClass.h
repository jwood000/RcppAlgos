#ifndef EXPOSE_CLASS_H
#define EXPOSE_CLASS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP StartOverGlue(SEXP ext);
    SEXP CombClassNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo,
                      SEXP Rparallel, SEXP RstdFun, SEXP Rrho,
                      SEXP R_RFunVal, SEXP RmainFun, SEXP RcompFun,
                      SEXP Rtarget, SEXP RKeepRes, SEXP Rtolerance,
                      SEXP RmIsNull, SEXP RretVal);
    SEXP NextCombGlue(SEXP ext);
    SEXP NextNumCombGlue(SEXP ext, SEXP Rnum);
    SEXP NextGatherGlue(SEXP ext);
    SEXP PrevCombGlue(SEXP ext);
    SEXP PrevNumCombGlue(SEXP ext, SEXP Rnum);
    SEXP PrevGatherGlue(SEXP ext);
    SEXP CurrCombGlue(SEXP ext);
    SEXP SourceVectorGlue(SEXP ext);
    SEXP RandomAccessGlue(SEXP ext, SEXP RIndexVec);
    SEXP FrontGlue(SEXP ext);
    SEXP BackGlue(SEXP ext);
    SEXP SummaryGlue(SEXP ext);
}

#endif
