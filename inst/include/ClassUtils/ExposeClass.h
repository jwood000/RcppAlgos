#ifndef EXPOSE_CLASS_H
#define EXPOSE_CLASS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    void StartOverGlue(SEXP ext);
    SEXP ComboNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo);
    SEXP ComboApplyNew(SEXP RVals, SEXP RboolVec, SEXP freqInfo,
                       SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal);
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
