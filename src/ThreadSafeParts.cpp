// #include "Combinations/NthCombination.h"
// #include "Combinations/ComboManager.h"
// #include <RcppThread/ThreadPool.hpp>
// #include <gmp.h>
// 
// template <typename T>
// void ThreadSafeGenPartitions(T* mat, const std::vector<T> &v, int n, int m,
//                              bool Parallel, bool IsRep, bool IsMult, bool IsGmp,
//                              const std::vector<int> &freqs, std::vector<int> z,
//                              const std::vector<int> &myReps, double lower,
//                              mpz_t lowerMpz, int nRows, int nThreads) {
// 
//     if (Parallel) {
//         RcppParallel::RMatrix<T> parMat(mat, nRows, m);
//         RcppThread::ThreadPool pool(nThreads);
//         
//         const int lastElem = lenV - 1;
//         const int lastCol = part.width - 1;
//         int boundary = lastCol;
//         int edge = boundary - 1;
// 
//         const int lastRow = nRows - 1;
//         const int stepSize = lastRow / nThreads;
//         int nextStep = stepSize;
//         int step = 0;
// 
//         const nthPartPtr nthPartFun = GetNthPartFunc(IsMult, IsRep, IsGmp);
//         std::vector<std::vector<int>> zs(nThreads, z);
// 
//         for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
//             pool.push(std::cref(PartsGenParallel<T>), std::ref(parMat), std::cref(v),
//                       std::ref(zs[j]), n, m, step, nextStep, std::cref(freqs),
//                       IsMult, IsRep);
// 
//             if (IsGmp) {
//                 mpz_add_ui(lowerMpz, lowerMpz, stepSize);
//             } else {
//                 lower += stepSize;
//             }
// 
//             zs[j + 1] = nthPartFun(n, m, lower, lowerMpz, myReps);
//         }
// 
//         pool.push(std::cref(PartsGenParallel<T>),
//                   std::ref(parMat), std::cref(v), std::ref(zs.back()),
//                   n, m, step, nRows, std::cref(freqs), IsMult, IsRep);
// 
//         pool.join();
//     } else {
//         PartsGenManager(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
//     }
// }
// 
// template void ThreadSafeGenPartitions(int*, const std::vector<int>&, int, int,
//                                       bool, bool, bool, bool, const std::vector<int>&,
//                                       std::vector<int>, const std::vector<int>&, double,
//                                       mpz_t, int, int);
// 
// template void ThreadSafeGenPartitions(double*, const std::vector<double>&, int, int,
//                                       bool, bool, bool, bool, const std::vector<int>&,
//                                       std::vector<int>, const std::vector<int>&, double,
//                                       mpz_t, int, int);
