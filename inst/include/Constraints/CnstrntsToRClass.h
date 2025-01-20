#pragma once

#include "Constraints/ConstraintsClass.h"
#include "ClassUtils/ComboResClass.h"

class CnstrntsToR : public ComboRes {
private:
    SEXP GetNext();
    SEXP GetNextN(int n);

    bool keepGoing = true;
    const int maxRows;

    int upperBoundInt;
    int upperBoundDbl;

    std::vector<int> currIntVec;
    std::vector<double> currDblVec;

    const std::vector<int> origTarIntVals;
    const std::vector<double> origTarVals;

    std::unique_ptr<ConstraintsClass<int>> CnstrtInt;
    std::unique_ptr<ConstraintsClass<double>> CnstrtDbl;

public:

    CnstrntsToR(
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
    SEXP nextNumIters(SEXP RNum);
    SEXP nextGather();
    SEXP currIter();
    SEXP summary();
};
