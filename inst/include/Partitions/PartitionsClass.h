#pragma once

#include "ClassUtils/ComboResClass.h"
#include "Partitions/NextPartition.h"
#include "Partitions/NthPartition.h"
#include "Sample/SamplePartitions.h"

class Partitions : public ComboRes {
private:
    int edge;
    int pivot;
    int tarDiff;
    int boundary;

    bool bAddOne;

    const bool paragon;
    const bool stdPartNext;
    const bool stdCompZeroSpesh;
    const bool genCompZeroSpesh;

    const int lastCol;
    const int lastElem;

    std::vector<int> rpsCnt;
    const nextPartsPtr nextParts;
    const nthPartsPtr nthParts;

    void SetPartValues();
    void MoveZToIndex();

    // We need a dedicated method for multisets to
    // avoid issues with updating the indexing vector
    SEXP MultisetMatrix(int nRows);

public:

    Partitions(
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
    SEXP randomAccess(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP summary();
};
