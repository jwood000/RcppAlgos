#include "Permutations/NthPermutation.h"
#include "Permutations/PermuteResGlue.h"
#include <RcppThread/ThreadPool.hpp>
#include <gmp.h>

template <typename T>
void PermResMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                 int n, int m, bool Parallel, bool IsRep, bool IsMult,
                 bool IsGmp, const std::vector<int> &freqs, std::vector<int> &z,
                 const std::vector<int> &myReps, double lower, mpz_t lowerMpz,
                 int nRows, int nThreads) {
    
    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, m);
        RcppThread::ThreadPool pool(nThreads);

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthPermPtr nthPermFun = GetNthPermFunc(IsMult, IsRep, IsGmp);
        std::vector<std::vector<int>> zs(nThreads, z);

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(PermuteResPar<T>), std::ref(parMat), std::cref(v),
                      std::ref(zs[j]), n, m, step, nextStep, std::cref(freqs),
                      std::cref(myFun), IsMult, IsRep);

            if (IsGmp) {
                mpz_add_ui(lowerMpz, lowerMpz, stepSize);
            } else {
                lower += stepSize;
            }

            SetStartPerm(zs[j + 1], nthPermFun, myReps, n,
                         m, lower, lowerMpz, IsRep, IsMult);
        }

        pool.push(std::cref(PermuteResPar<T>),
                  std::ref(parMat), std::cref(v), std::ref(zs.back()), n, m,
                  step, nRows, std::cref(freqs), std::cref(myFun), IsMult, IsRep);

        pool.join();
    } else {
        PermuteResStd(mat, v, z, n, m, nRows, IsMult, IsRep, freqs, myFun);
    }
}

template void PermResMain(int*, const std::vector<int>&,
                          const funcPtr<int>, int, int, bool, bool, bool,
                          bool, const std::vector<int>&, std::vector<int>&,
                          const std::vector<int>&, double, mpz_t, int, int);

template void PermResMain(double*, const std::vector<double>&,
                          const funcPtr<double>, int, int, bool, bool, bool,
                          bool, const std::vector<int>&, std::vector<int>&,
                          const std::vector<int>&, double, mpz_t, int, int);
