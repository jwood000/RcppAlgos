#ifndef COMBO_CLASS_H
#define COMBO_CLASS_H

#include "ClassUtils/NextCombinatorics.h"
#include "Permutations/NthPermutation.h"
#include "ClassUtils/GetPrevCombPerm.h"
#include "Sample/SampCombPermStd.h"
#include "ClassUtils/ClassUtils.h"
#include "Cpp14MakeUnique.h"
#include "ImportExportMPZ.h"
#include "GetCombPerm.h"
#include "SetUpUtils.h"
#include "NthResult.h"

class Combo {
private:
    SEXP VecReturn();
    SEXP MatForward(int nRows);
    SEXP MatReverse(int nRows);

protected:
    const int n;
    const int m;
    const int m1;
    int RTYPE;
    const int maxThreads;

    const SEXP sexpVec;
    const SEXP sexpNThreads;

    // IsGmp may change depending on whether we have partitions
    bool IsGmp;
    const bool IsFactor;
    const bool IsComb;
    const bool IsMult;
    const bool IsRep;
    const bool Parallel;

    double computedRows;
    const VecType myType;
    mpz_t computedRowsMpz;

    double dblTemp;
    mpz_t mpzTemp;

    std::vector<int> z;
    std::vector<int> vInt;
    std::vector<double> vNum;

    const std::vector<int> freqs;
    std::vector<int> myReps;

    // This has to be initialized later becuase it
    // depends on freqs.size, IsMult, and n
    const int n1;

    double dblIndex;
    mpz_t mpzIndex;

    SEXP myClass;
    SEXP myLevels;

    const nthResultPtr nthResFun;
    const nextIterPtr nextIter;
    const prevIterPtr prevIter;

    bool prevIterAvailable;
    SEXP ToSeeLast(bool AdjustIdx = true);
    SEXP ToSeeFirst(bool AdjustIdx = true);

public:

    Combo(
        SEXP Rv, int Rm, SEXP RcompRow, const std::vector<int> &bVec,
        const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
        const std::vector<int> &RvInt, const std::vector<double> &RvNum,
        VecType typePass, int RmaxThreads, SEXP RnThreads, bool Rparallel
    );

    virtual ~Combo() = default;
    virtual void startOver();
    virtual SEXP nextComb();
    virtual SEXP prevComb();
    virtual SEXP nextNumCombs(SEXP RNum);
    virtual SEXP prevNumCombs(SEXP RNum);
    virtual SEXP nextGather();
    virtual SEXP prevGather();
    virtual SEXP currComb();
    virtual SEXP randomAccess(SEXP RindexVec);
    virtual SEXP front();
    virtual SEXP back();
    SEXP sourceVector() const;
    virtual SEXP summary();
};

#endif
