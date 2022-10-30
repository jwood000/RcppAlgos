#pragma once

#include "ClassUtils/GetPrevCombPermApply.h"
#include "ClassUtils/ComboClass.h"
#include "Sample/SampleApply.h"
#include "GetCombPermApply.h"

class ComboApply : public Combo {
private:
    const SEXP rho;
    const SEXP stdFun;
    const SEXP RFunVal;

    SEXP VecApplyReturn();
    SEXP ApplyForward(int nRows);
    SEXP ApplyReverse(int nRows);

public:
    ComboApply(
        SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
        SEXP RstdFun, SEXP Rrho, SEXP R_RFunVal
    );

    void startOver();
    SEXP nextComb();
    SEXP prevComb();
    SEXP nextNumCombs(SEXP RNum);
    SEXP prevNumCombs(SEXP RNum);
    SEXP nextGather();
    SEXP prevGather();
    SEXP currComb();
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
};
