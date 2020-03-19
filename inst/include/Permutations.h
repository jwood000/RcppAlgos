#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "NextStandard.h"
#include "StandardCount.h"
#include "Cpp14MakeUnique.h"
#include <RcppThread/ThreadPool.hpp>
#include <cmath>

constexpr int unrollSize = 8;

template <typename typeMatrix, typename typeVector>
void RepUnroller(typeMatrix &matRcpp, const typeVector &v, 
                 int strt, int last, int ind, int lastUnroll) {
    
    for (int i = strt; i < lastUnroll; i += unrollSize) {
        matRcpp(i, 0) = v[ind];
        matRcpp(i + 1, 0) = v[ind];
        matRcpp(i + 2, 0) = v[ind];
        matRcpp(i + 3, 0) = v[ind];
        matRcpp(i + 4, 0) = v[ind];
        matRcpp(i + 5, 0) = v[ind];
        matRcpp(i + 6, 0) = v[ind];
        matRcpp(i + 7, 0) = v[ind];
    }
    
    for (int i = lastUnroll; i < last; ++i)
        matRcpp(i, 0) = v[ind];
}

template <typename typeMatrix, typename typeVector>
void StandardUnroller(typeMatrix &matRcpp, const int *const indexMat, const typeVector &v,
                      int m, int strt, int last, int first, int lastUnroll) {
    
    for (int j = first, k = 0; j < m; ++j) {
        for (int i = strt; i < lastUnroll; i += unrollSize, k += unrollSize) {
            matRcpp(i, j) = v[indexMat[k]];
            matRcpp(i + 1, j) = v[indexMat[k + 1]];
            matRcpp(i + 2, j) = v[indexMat[k + 2]];
            matRcpp(i + 3, j) = v[indexMat[k + 3]];
            matRcpp(i + 4, j) = v[indexMat[k + 4]];
            matRcpp(i + 5, j) = v[indexMat[k + 5]];
            matRcpp(i + 6, j) = v[indexMat[k + 6]];
            matRcpp(i + 7, j) = v[indexMat[k + 7]];
        }
        
        for (int i = lastUnroll; i < last; ++i, ++k)
            matRcpp(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteWorker(typeMatrix &matRcpp, const int *const indexMat, typeVector v,
                   int m, int strt, int last, int ind, int first, int unrollRem, bool IsRep) {
    
    const int lastUnroll = last - unrollRem;
    // For IsRep case, we are not setting the first column because we know that it
    // simply increments. This is taken into account in the indexMat preparation.
    if (IsRep) {RepUnroller(matRcpp, v, strt, last, ind, lastUnroll);}
    StandardUnroller(matRcpp, indexMat, v, m, strt, last, first, lastUnroll);
}

template <typename typeMatrix, typename typeVector>
void PermuteLoadIndex(typeMatrix &matRcpp, int *const indexMat, const typeVector &v,
                      std::vector<int> &z, int n, int m, int segment, bool IsRep) {
    
    const int maxInd = n - 1;
    const int lastCol = m - 1;
    
    if (IsRep) {
        for (int count = 0; count < segment; ++count) {
            // N.B. In PermuteGeneral we start j at 0
            for (int j = 1, k = count; j < m; ++j, k += segment) {
                matRcpp(count, j) = v[z[j]];
                indexMat[k] = z[j];
            }
            
            matRcpp(count, 0) = v[z.front()];
            
            // N.B. In PermuteGeneral we decrement i until i == 0
            for (int i = lastCol; i > 0; --i) {
                if (z[i] != maxInd) {
                    ++z[i];
                    break;
                } else {
                    z[i] = 0;
                }
            }
        }
    } else {
        auto arrPerm = FromCpp14::make_unique<int[]>(n);
        
        for (int i = 0; i < n; ++i)
            arrPerm[i] = z[i];
        
        if (m == n) {
            for (int count = 0; count < segment; ++count) {
                for (int j = 0, k = count; j < m; ++j, k += segment) {
                    matRcpp(count, j) = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }
                
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (int count = 0; count < segment; ++count) {
                for (int j = 0, k = count; j < m; ++j, k += segment) {
                    matRcpp(count, j) = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }
                
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteParallel(typeMatrix &matRcpp, typeVector v, std::vector<int> z, 
                     int n, int m, int nRows, int segment, int nThreads, bool IsRep) {

    const int first = (IsRep) ? 1 : 0;
    const std::size_t indexMatSize = static_cast<std::size_t>(segment)
                                     * static_cast<std::size_t>(m - first);
    
    auto indexMat = FromCpp14::make_unique<int[]>(indexMatSize);
    PermuteLoadIndex(matRcpp, indexMat.get(), v, z, n, m, segment, IsRep);
    
    int ind = 1;
    int strt = segment;
    int last = strt + segment;
    int unrollRem = segment % unrollSize;
    RcppThread::ThreadPool pool(nThreads);

    for (; last <= nRows;) {
        for (int j = 0; j < nThreads && last <= nRows; ++j, strt += segment, last += segment, ++ind) {
            if (!IsRep) {std::swap(v[0], v[ind]);}
            pool.push(std::cref(PermuteWorker<typeMatrix, typeVector>), std::ref(matRcpp),
                      indexMat.get(), v, m, strt, last, ind, first, unrollRem, IsRep);
        }
        
        if (last <= nRows) {pool.wait();}
    }

    pool.join();

    if (ind < n && strt < nRows) {
        const int skip = last - nRows;
        unrollRem = nRows % unrollSize;
        const int lastUnroll = nRows - unrollRem;

        if (IsRep) {
            RepUnroller(matRcpp, v, strt, nRows, ind, lastUnroll);
        } else {
            std::swap(v[0], v[ind]);
        }

        for (int j = first, k = 0; j < m; ++j, k += skip)
            for (int i = strt; i < nRows; ++i, ++k)
                matRcpp(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteSerialNoRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                        int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    constexpr int first = 0;
    const int segment = NumPermsNoRep(n - 1, m - 1);
    const std::size_t indexMatSize = static_cast<std::size_t>(segment)
                                     * static_cast<std::size_t>(m - first);
    
    auto indexMat = FromCpp14::make_unique<int[]>(indexMatSize);
    PermuteLoadIndex(matRcpp, indexMat.get(), v, z, n, m, segment, false);

    int ind = 1;
    strt = segment;
    int last = strt + segment;
    int unrollRem = segment % unrollSize;
    typeVector vCopy(v);

    for (; last <= nRows; strt += segment, last += segment, ++ind) {
        typeVector vTemp(1);
        vTemp[0] = vCopy[0];
        vCopy[0] = vCopy[ind];
        vCopy[ind] = vTemp[0];

        PermuteWorker(matRcpp, indexMat.get(), vCopy, m,
                      strt, last, ind, first, unrollRem, false);
    }

    if (ind < static_cast<int>(v.size()) && strt < nRows) {
        const int skip = last - nRows;
        unrollRem = nRows % unrollSize;

        typeVector vTemp(1);
        vTemp[0] = vCopy[0];
        vCopy[0] = vCopy[ind];
        vCopy[ind] = vTemp[0];

        for (int j = first, k = 0; j < m; ++j, k += skip)
            for (int i = strt; i < nRows; ++i, ++k)
                matRcpp(i, j) = vCopy[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteSerialRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                      int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    constexpr int first = 1;
    const int segment = std::pow(static_cast<double>(n),
                                 static_cast<double>(m - 1));
    const std::size_t indexMatSize = static_cast<std::size_t>(segment)
                                     * static_cast<std::size_t>(m - first);
    
    auto indexMat = FromCpp14::make_unique<int[]>(indexMatSize);
    PermuteLoadIndex(matRcpp, indexMat.get(), v, z, n, m, segment, true);
    
    int ind = 1;
    strt = segment;
    int last = strt + segment;
    int unrollRem = segment % unrollSize;
    
    for (; last <= nRows; strt += segment, last += segment, ++ind) {
        PermuteWorker(matRcpp, indexMat.get(), v, m, 
                      strt, last, ind, first, unrollRem, true);
    }
    
    if (ind < static_cast<int>(v.size()) && strt < nRows) {
        const int skip = last - nRows;
        unrollRem = nRows % unrollSize;
        const int lastUnroll = nRows - unrollRem;
        RepUnroller(matRcpp, v, strt, nRows, ind, lastUnroll);
        
        for (int j = first, k = 0; j < m; ++j, k += skip)
            for (int i = strt; i < nRows; ++i, ++k)
                matRcpp(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteGeneralNoRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                         int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    const int maxInd = n - 1;
    const int numR1 = nRows - 1;
    auto arrPerm = FromCpp14::make_unique<int[]>(n);
    
    for (int i = 0; i < n; ++i)
        arrPerm[i] = z[i];
    
    if (m == n) {
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        const int lastCol = m - 1;
        
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j)
        matRcpp(numR1, j) = v[arrPerm[j]];
}

template <typename typeMatrix, typename typeVector>
void PermuteGeneralRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                       int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    const int maxInd = n - 1;
    const int lastCol = m - 1;
    
    for (int count = strt; count < nRows; ++count) {
        for (int j = 0; j < m; ++j)
            matRcpp(count, j) = v[z[j]];
        
        for (int i = lastCol; i >= 0; --i) {
            if (z[i] != maxInd) {
                ++z[i];
                break;
            } else {
                z[i] = 0;
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermutation(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                         int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    const int lenFreqs = z.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);
    
    const int numR1 = nRows - 1;
    const int maxInd = lenFreqs - 1;
    
    for (int j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (m == lenFreqs) {
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];
            
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        const int lastCol = m - 1;
        
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];
            
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }
    
    // Get last permutation
    for (int j = 0; j < m; ++j)
        matRcpp(numR1, j) = v[arrPerm[j]];
}

#endif
