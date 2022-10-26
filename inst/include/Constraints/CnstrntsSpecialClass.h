#pragma once

#include "ClassUtils/ComboResClass.h"

class CnstrntsSpecial : public ComboRes {
private:
    int count;
    bool keepGoing;

public:
    CnstrntsSpecial(
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
    SEXP nextComb();
    SEXP nextNumCombs(SEXP RNum);
    SEXP nextGather();
    SEXP currComb();
    SEXP summary();
};
