#include "Permutations/PermuteManager.h"
#include "NthResult.h"
#include "ThreadBatch.h"
#include <gmpxx.h>

template <typename T>
void ThreadSafePermutations(T* mat, const std::vector<T> &v, int n, int m,
                            int phaseOne, bool generalRet, bool Parallel,
                            bool IsRep, bool IsMult, bool IsGmp,
                            const std::vector<int> &freqs, std::vector<int> &z,
                            const std::vector<int> &myReps, double lower,
                            mpz_class &lowerMpz, int nRows,
                            int nThreads) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, m);

        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthResultPtr nthResFun = GetNthResultFunc(false, IsMult,
                                                        IsRep, IsGmp);
        std::vector<std::vector<int>> zs(nThreads, z);
        RcppAlgos::ThreadBatch batch(nThreads);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            batch.launch([&, j, step, nextStep] {
                PermuteParallel(parMat, v, zs[j], n, m, step, nextStep,
                                freqs, IsMult, IsRep);
            });

            SetNextIter(myReps, zs[j + 1], nthResFun, lower, lowerMpz,
                        stepSize, n, m, IsGmp, false, IsRep, IsMult);
        }

        batch.launch([&, step] {
            PermuteParallel(parMat, v, zs.back(), n, m, step, nRows,
                            freqs, IsMult, IsRep);
        });

        batch.finish();
    } else {
        PermuteManager(mat, v, z, n, m, nRows, phaseOne,
                       generalRet, IsMult, IsRep, freqs);
    }
}

template void ThreadSafePermutations(
    int*, const std::vector<int>&, int, int, int, bool, bool, bool, bool,
    bool, const std::vector<int>&, std::vector<int>&, const std::vector<int>&,
    double, mpz_class&, int, int
);

template void ThreadSafePermutations(
    double*, const std::vector<double>&, int, int, int, bool, bool, bool,
    bool, bool, const std::vector<int>&, std::vector<int>&,
    const std::vector<int>&, double, mpz_class&, int, int
);
