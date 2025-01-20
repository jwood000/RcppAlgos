#pragma once

#include "Constraints/GetContraints.h"
#include "ClassUtils/ComboClass.h"

class ComboRes : public Combo {
protected:
    SEXP ApplyFun(SEXP res);
    SEXP VecReturn();
    SEXP MatrixReturn(int nRows);

    const int cap;
    const int width;
    const int nCols;
    const int strtLen;

    bool bLower;
    bool bUpper;

    const bool KeepRes;
    const bool numUnknown;

    const double cnstrtCount;
    const mpz_class cnstrtCountMpz;

    std::vector<int> tarIntVals;
    std::vector<double> tarVals;

    const ConstraintType ctype;
    const PartDesign part;

    const std::string mainFun;
    const std::string funTest;

    const std::vector<std::string> compVec;
    const funcPtr<double> funDbl;
    const funcPtr<int> funInt;

public:

    ComboRes(
        SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
        const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
        std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
        std::vector<int> &RstartZ, const std::string &RmainFun,
        const std::string &RFunTest, funcPtr<double> RfunDbl,
        ConstraintType Rctype, int RstrtLen, int Rcap, bool RKeepRes,
        bool RnumUnknown, double RcnstrtRows, const mpz_class &RcnstrtRowsMpz
    );

    void startOver();
    SEXP nextIter();
    SEXP prevIter();
    SEXP nextNumIters(SEXP RNum);
    SEXP prevNumIters(SEXP RNum);
    SEXP nextGather();
    SEXP prevGather();
    SEXP currIter();
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP summary();
};
