#ifndef CONSTRAINTS_CLASS_H
#define CONSTRAINTS_CLASS_H

#include "Constraints/ConstraintsUtils.h"
#include "Constraints/NextGeneralRes.h"
#include "ClassUtils/ComboResClass.h"

template <typename T>
class Constraints : public ComboRes {
private:
    void PrepareConstraints();

    const int m2;
    const int nMinusM;
    const int maxZ;
    const int freqsSize;
    const int pentExtreme;

    bool check_0;
    bool check_1;

    const funcPtr<T> fun;
    const compPtr<T> comp;
    std::vector<int> zIndex;

    std::vector<T> v;
    std::vector<T> testVec;

    const std::vector<T> targetVals;
    const nextCnstrtPtr<T> nextCnstrnt;

public:

    Constraints(
        SEXP Rv, int Rm, SEXP RcompRows, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnumThreads, bool Rparallel,
        const PartDesign &Rpart, const std::vector<std::string> &RcompVec,
        std::vector<double> &RtarVals, std::vector<int> &RtarIntVals,
        std::vector<int> &RstartZ, const std::string &RmainFun,
        funcPtr<double> RfunDbl, ConstraintType Rctype, int RstrtLen,
        int Rcap, bool RKeepRes, bool RnumUnknown, double RcnstrtRows,
        mpz_t RcnstrtRowsMpz
    );

    void startOver();
    SEXP nextComb();
    SEXP nextNumCombs(SEXP RNum);
    SEXP nextGather();
    SEXP currComb();
    SEXP summary();
};

#endif
