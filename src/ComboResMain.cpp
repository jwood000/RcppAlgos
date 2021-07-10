#include "Combinations/ComboResGlue.h"
#include "NthResult.h"
#include <thread>
#include <gmp.h>

template <typename T>
void ComboResMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                  int n, int m, bool Parallel, bool IsRep, bool IsMult,
                  bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_t lowerMpz, int nRows, int nThreads) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, m);
        std::vector<std::thread> threads;

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthResultPtr nthResFun = GetNthResultFunc(true, IsMult,
                                                        IsRep, IsGmp);
        std::vector<std::vector<int>> zs(nThreads, z);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            threads.emplace_back(std::cref(ComboResPar<T>), std::ref(parMat),
                                 std::cref(v), std::ref(zs[j]), n, m, step,
                                 nextStep, std::cref(freqs), std::cref(myFun),
                                 IsMult, IsRep);

            SetNextIter(myReps, zs[j + 1], nthResFun, lower, lowerMpz,
                        stepSize, n, m, IsGmp, true, IsRep, IsMult);
        }

        threads.emplace_back(std::cref(ComboResPar<T>), std::ref(parMat),
                             std::cref(v), std::ref(zs.back()), n, m, step,
                             nRows, std::cref(freqs), std::cref(myFun),
                             IsMult, IsRep);

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        ComboResStd(mat, v, z, n, m, nRows, IsMult, IsRep, freqs, myFun);
    }
}

template void ComboResMain(int*, const std::vector<int>&,
                           const funcPtr<int>, int, int, bool, bool, bool,
                           bool, const std::vector<int>&, std::vector<int>&,
                           const std::vector<int>&, double, mpz_t, int, int);

template void ComboResMain(double*, const std::vector<double>&,
                           const funcPtr<double>, int, int, bool, bool, bool,
                           bool, const std::vector<int>&, std::vector<int>&,
                           const std::vector<int>&, double, mpz_t, int, int);
