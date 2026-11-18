#include "Combinations/ComboManager.h"
#include "NthResult.h"
#include "ThreadBatch.h"
#include <gmpxx.h>

template <typename T>
void ThreadSafeCombinations(T* mat, const std::vector<T> &v, int n, int m,
                            bool Parallel, bool IsRep, bool IsMult, bool IsGmp,
                            const std::vector<int> &freqs, std::vector<int> &z,
                            const std::vector<int> &myReps, double lower,
                            mpz_class &lowerMpz, int nRows, int nThreads) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, m);

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthResultPtr nthResFun = GetNthResultFunc(true, IsMult,
                                                        IsRep, IsGmp);
        std::vector<std::vector<int>> zs(nThreads, z);
        RcppAlgos::ThreadBatch batch(nThreads);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            batch.launch([&, j, step, nextStep] {
                ComboParallel(parMat, v, zs[j], n, m, step, nextStep,
                              freqs, IsMult, IsRep);
            });

            SetNextIter(myReps, zs[j + 1], nthResFun, lower, lowerMpz,
                        stepSize, n, m, IsGmp, true, IsRep, IsMult);
        }

        batch.launch([&, step] {
            ComboParallel(parMat, v, zs.back(), n, m, step, nRows,
                          freqs, IsMult, IsRep);
        });

        batch.finish();
    } else {
        ComboManager(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    }
}

template void ThreadSafeCombinations(int*, const std::vector<int>&, int, int,
                                     bool, bool, bool, bool, const std::vector<int>&,
                                     std::vector<int>&, const std::vector<int>&, double,
                                     mpz_class&, int, int);

template void ThreadSafeCombinations(double*, const std::vector<double>&, int, int,
                                     bool, bool, bool, bool, const std::vector<int>&,
                                     std::vector<int>&, const std::vector<int>&, double,
                                     mpz_class&, int, int);
