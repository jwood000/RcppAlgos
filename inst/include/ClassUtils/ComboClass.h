#ifndef COMBO_CLASS_H
#define COMBO_CLASS_H

#include "ClassUtils/NextCombinatorics.h"
#include "Permutations/NthPermutation.h"
#include "ClassUtils/ComboClassUtils.h"
#include "ClassUtils/GetPrevCombPerm.h"
#include "Sample/SampCombPermStd.h"
#include "Cpp14MakeUnique.h"
#include "ImportExportMPZ.h"
#include "GetCombPerm.h"
#include "SetUpUtils.h"
#include "NthResult.h"

class Combo {
private:
    const int n;
    const int m;
    const int m1;
    const int RTYPE;
    const int maxThreads;
    const SEXP sexpVec;
    const bool IsFactor;
    const bool IsComb;
    const bool IsMult;
    const bool IsRep;
    const bool IsGmp;
    const bool Parallel;
    double computedRows;
    const VecType myType;
    mpz_t computedRowsMpz[1];
    
    double dblTemp;
    mpz_t mpzTemp;
    std::vector<int> z;
    const std::vector<int> vInt;
    const std::vector<double> vNum;
    
    const std::vector<int> freqs;
    const std::vector<int> myReps;
    
    double dblIndex;
    mpz_t mpzIndex[1];

    SEXP myClass;
    SEXP myLevels;
    SEXP SexpNThreads;
    
    const nthResultPtr nthResFun;
    const nextIterPtr nextIter;
    const prevIterPtr prevIter;
    
    // This has to be initialized later becuase it
    // depends on freqs.size, IsMult, and n
    const int n1;
    
    SEXP VecReturn();
    SEXP MatForward(int nRows);
    SEXP MatReverse(int nRows);

public:
    
    Combo(SEXP Rv, int Rm, SEXP RcompRow, const std::vector<int> &bVec,
          const std::vector<int> &Rreps, const std::vector<int> &Rfreqs,
          const std::vector<int> &RvInt, const std::vector<double> &RvNum,
          VecType typePass, int RmaxThreads, SEXP RnumThreads);

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
    SEXP sourceVector() const;
    SEXP summary();
};

#endif
