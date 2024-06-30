#pragma once

#include "ClassUtils/NextCombinatorics.h"
#include "Permutations/NthPermutation.h"
#include "ClassUtils/GetPrevCombPerm.h"
#include "Sample/SampCombPermStd.h"
#include "ClassUtils/ClassUtils.h"
#include "ClassUtils/Iterator.h"
#include "GetCombPerm.h"
#include "NthResult.h"

class Combo : public Iterator {
private:

    SEXP MatForward(int nRows, int numIncrement);
    SEXP MatReverse(int nRows);

protected:

    const int m;
    const int m1;

    const bool IsFactor;
    const bool IsComb;
    const bool IsMult;
    const bool IsRep;
    const VecType myType;

    double dblTemp;
    mpz_class mpzTemp;

    std::vector<int> z;
    std::vector<int> vInt;
    std::vector<double> vNum;

    const std::vector<int> freqs;
    std::vector<int> myReps;

    // This has to be initialized later becuase it
    // depends on freqs.size, IsMult, and n
    const int n1;

    double dblIndex;
    mpz_class mpzIndex;

    SEXP myClass;
    SEXP myLevels;

    const nthResultPtr nthResFun;
    const nextIterPtr nextIter;
    const prevIterPtr prevIter;

    bool prevIterAvailable;

    SEXP BasicVecReturn();
    SEXP ToSeeLast(bool AdjustIdx = true);
    SEXP ToSeeFirst(bool AdjustIdx = true);

public:

    Combo(
        SEXP Rv, int Rm, SEXP RcompRow, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnThreads, bool Rparallel
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
    SEXP summary();
};
