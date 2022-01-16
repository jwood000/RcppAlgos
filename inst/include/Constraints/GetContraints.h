#ifndef GET_CONSTRAINTS_H
#define GET_CONSTRAINTS_H

#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsTypes.h"

SEXP GetConstraints(
    const PartDesign &part, const std::vector<std::string> &compVec,
    const std::vector<int> &freqs, std::vector<int> &myReps,
    std::vector<double> &vNum, std::vector<int> &vInt,
    std::vector<double> &tarVals, std::vector<int> &tarIntVals,
    std::vector<int> &startZ, const std::string &mainFun,
    const std::string &funTest, funcPtr<double> funDbl, double lower,
    mpz_t lowerMpz, double userNum, ConstraintType ctype, VecType myType,
    int nThreads, int nRows, int n, int strtLen, int cap, int m,
    bool IsComb, bool Parallel, bool IsGmp, bool IsRep, bool IsMult,
    bool bUpper, bool KeepRes, bool numUnknown
);

#endif
