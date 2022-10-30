#include "Permutations/ThreadSafePerm.h"
#include "Combinations/ThreadSafeComb.h"
#include "Permutations/PermuteManager.h"
#include "Combinations/ComboManager.h"
#include "SetUpUtils.h"

void CharacterGlue(SEXP mat, SEXP v, bool IsComb,
                   std::vector<int> &z, int n, int m, int nRows,
                   const std::vector<int> &freqs, bool IsMult, bool IsRep) {

    if (IsComb) {
        ComboCharacter(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    } else {
        PermuteCharacter(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    }
}

template <typename T>
void ManagerGlue(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int n, int m, int nRows, bool IsComb, int phaseOne,
                 bool generalRet, const std::vector<int> &freqs,
                 bool IsMult, bool IsRep) {

    if (IsComb) {
        ComboManager(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    } else {
        PermuteManager(mat, v, z, n, m, nRows, phaseOne,
                       generalRet, IsMult, IsRep, freqs);
    }
}

template <typename T>
void ParallelGlue(T* mat, const std::vector<T> &v, int n, int m, int phaseOne,
                  bool generalRet, bool IsComb, bool Parallel, bool IsRep,
                  bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_class lowerMpz, int nRows, int nThreads) {

    if (IsComb) {
        ThreadSafeCombinations(mat, v, n, m, Parallel, IsRep,
                               IsMult, IsGmp, freqs, z, myReps,
                               lower, lowerMpz, nRows, nThreads);
    } else {
        ThreadSafePermutations(mat, v, n, m, phaseOne, generalRet, Parallel,
                               IsRep, IsMult, IsGmp, freqs, z, myReps, lower,
                               lowerMpz, nRows, nThreads);
    }
}

SEXP GetCombPerms(SEXP Rv, const std::vector<double> &vNum,
                  const std::vector<int> &vInt, int n, int m, int phaseOne,
                  bool generalRet, bool IsComb, bool Parallel, bool IsRep,
                  bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_class &lowerMpz, int nRows,
                  int nThreads, VecType myType) {

    switch (myType) {
        case VecType::Character : {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp res = Rf_allocMatrix(STRSXP, nRows, m);

            CharacterGlue(res, charVec, IsComb, z, n, m,
                          nRows, freqs, IsMult, IsRep);

            return res;
        } case VecType::Complex : {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i) {
                stlCmplxVec[i] = vecCmplx[i];
            }

            cpp11::sexp res = Rf_allocMatrix(CPLXSXP, nRows, m);
            Rcomplex* matCmplx = COMPLEX(res);

            ManagerGlue(matCmplx, stlCmplxVec, z, n, m, nRows, IsComb,
                        phaseOne, generalRet, freqs, IsMult, IsRep);

            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* rawVec = RAW(Rv);

            for (int i = 0; i < n; ++i) {
                stlRawVec[i] = rawVec[i];
            }

            cpp11::sexp res = Rf_allocMatrix(RAWSXP, nRows, m);
            Rbyte* rawMat = RAW(res);

            ManagerGlue(rawMat, stlRawVec, z, n, m, nRows, IsComb,
                        phaseOne, generalRet, freqs, IsMult, IsRep);

            return res;
        } case VecType::Logical : {
            std::vector<int> vBool(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vBool[i] = vecBool[i];
            }

            cpp11::sexp res = Rf_allocMatrix(LGLSXP, nRows, m);
            int* matBool = LOGICAL(res);

            ManagerGlue(matBool, vBool, z, n, m, nRows, IsComb,
                        phaseOne, generalRet, freqs, IsMult, IsRep);

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, nRows, m);
            int* matInt = INTEGER(res);

            ParallelGlue(matInt, vInt, n, m, phaseOne, generalRet, IsComb,
                         Parallel, IsRep, IsMult, IsGmp, freqs, z, myReps,
                         lower, lowerMpz, nRows, nThreads);

            if (Rf_isFactor(Rv)) SetFactorClass(res, Rv);
            return res;
        } default : {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, m);
            double* matNum = REAL(res);

            ParallelGlue(matNum, vNum, n, m, phaseOne, generalRet, IsComb,
                         Parallel, IsRep, IsMult, IsGmp, freqs, z, myReps,
                         lower, lowerMpz, nRows, nThreads);

            return res;
        }
    }
}
