#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "CombPermUtils.h"
#include "Cpp14MakeUnique.h"
#include <RcppThread.h>

constexpr std::size_t unrollSize = 8u;

template <typename typeMatrix, typename typeVector>
void RepUnroller(std::size_t r, const typeVector &v, const int *const indexMat,
                 typeMatrix &permuteMatrix, std::size_t start, std::size_t lastUnroll,
                 std::size_t last, std::size_t ind) {
    
    for (std::size_t i = start; i < lastUnroll; i += unrollSize) {
        permuteMatrix(i, 0u) = v[ind];
        permuteMatrix(i + 1u, 0u) = v[ind];
        permuteMatrix(i + 2u, 0u) = v[ind];
        permuteMatrix(i + 3u, 0u) = v[ind];
        permuteMatrix(i + 4u, 0u) = v[ind];
        permuteMatrix(i + 5u, 0u) = v[ind];
        permuteMatrix(i + 6u, 0u) = v[ind];
        permuteMatrix(i + 7u, 0u) = v[ind];
    }
    
    for (std::size_t i = lastUnroll; i < last; ++i)
        permuteMatrix(i, 0u) = v[ind];
}

template <typename typeMatrix, typename typeVector>
void StandardUnroller(std::size_t r, const typeVector &v, const int *const indexMat,
                      typeMatrix &permuteMatrix, std::size_t start, std::size_t lastUnroll,
                      std::size_t last, std::size_t ind, std::size_t first) {
    
    for (std::size_t j = first, k = 0u; j < r; ++j) {
        for (std::size_t i = start; i < lastUnroll; i += unrollSize, k += unrollSize) {
            permuteMatrix(i, j) = v[indexMat[k]];
            permuteMatrix(i + 1u, j) = v[indexMat[k + 1u]];
            permuteMatrix(i + 2u, j) = v[indexMat[k + 2u]];
            permuteMatrix(i + 3u, j) = v[indexMat[k + 3u]];
            permuteMatrix(i + 4u, j) = v[indexMat[k + 4u]];
            permuteMatrix(i + 5u, j) = v[indexMat[k + 5u]];
            permuteMatrix(i + 6u, j) = v[indexMat[k + 6u]];
            permuteMatrix(i + 7u, j) = v[indexMat[k + 7u]];
        }
        
        for (std::size_t i = lastUnroll; i < last; ++i, ++k)
            permuteMatrix(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteWorker(std::size_t r, bool IsRep, const typeVector &v, const int *const indexMat,
                   typeMatrix &permuteMatrix, std::size_t start, std::size_t last, 
                   std::size_t ind, std::size_t first, std::size_t unrollRem) {
    
    const std::size_t lastUnroll = last - unrollRem;
    
    // For IsRep case, we are not setting the first column because we know that it
    // simply increments. This is taken into account in the indexMat preparation.
    if (IsRep)
        RepUnroller(r, v, indexMat, permuteMatrix, start, lastUnroll, last, ind);
    
    StandardUnroller(r, v, indexMat, permuteMatrix, start, lastUnroll, last, ind, first);
}

template <typename typeMatrix, typename typeVector>
void PermuteLoadIndex(std::size_t n, std::size_t r, const typeVector &v, 
                      bool IsRep, std::size_t segment, std::vector<int> &z,
                      typeMatrix &permuteMatrix, int *const indexMat) {
    
    const std::size_t maxInd = n - 1u;
    const std::size_t lastCol = r - 1u;
    
    if (IsRep) {
        const int maxIndInt = maxInd;
        
        for (std::size_t count = 0u; count < segment; ++count) {
            // N.B. In PermuteGeneral we start j at 0u
            for (std::size_t j = 1u, k = count; j < r; ++j, k += segment) {
                permuteMatrix(count, j) = v[z[j]];
                indexMat[k] = z[j];
            }
            
            permuteMatrix(count, 0u) = v[z[0u]];
            
            // N.B. In PermuteGeneral we decrement k until k == 0
            for (int k = lastCol; k > 0; --k) {
                if (z[k] != maxIndInt) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }
    } else {
        auto arrPerm = FromCpp14::make_unique<int[]>(n);
        
        for (std::size_t i = 0u; i < n; ++i)
            arrPerm[i] = z[i];
        
        if (r == n) {
            for (std::size_t count = 0u; count < segment; ++count) {
                for (std::size_t j = 0u, k = count; j < r; ++j, k += segment) {
                    permuteMatrix(count, j) = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }
                
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (std::size_t count = 0u; count < segment; ++count) {
                for (std::size_t j = 0u, k = count; j < r; ++j, k += segment) {
                    permuteMatrix(count, j) = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }
                
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteParallel(std::size_t n, std::size_t r, typeVector v, bool IsRep, 
                     std::size_t uRowN, std::size_t segment, std::vector<int> &z,
                     typeMatrix &permuteMatrix, int nThreads) {

    const std::size_t first = (IsRep) ? 1u : 0u;
    auto indexMat = FromCpp14::make_unique<int[]>(segment * (r - first));
    PermuteLoadIndex(n, r, v, IsRep, segment, z, permuteMatrix, indexMat.get());
    
    std::size_t unrollRem = segment % unrollSize;
    std::size_t ind = 1u;
    std::size_t start = segment;
    std::size_t last = start + segment;
    RcppThread::ThreadPool pool(nThreads);

    for (; last <= uRowN;) {
        for (int j = 0; j < nThreads && last <= uRowN; ++j, start += segment, last += segment, ++ind) {
            if (!IsRep) {std::swap(v[0u], v[ind]);}
            pool.push(std::cref(PermuteWorker<typeMatrix, typeVector>), r, IsRep, v,
                      indexMat.get(), std::ref(permuteMatrix), start, last, ind, first, unrollRem);
        }
        if (last <= uRowN) {pool.wait();}
    }

    pool.join();

    if (ind < n && start < uRowN) {
        const std::size_t skip = last - uRowN;
        unrollRem = uRowN % unrollSize;
        const std::size_t lastUnroll = uRowN - unrollRem;

        if (IsRep) {
            RepUnroller(r, v, indexMat.get(), permuteMatrix, start, lastUnroll, uRowN, ind);
        } else {
            std::swap(v[0u], v[ind]);
        }

        for (std::size_t j = first, k = 0; j < r; ++j, k += skip)
            for (std::size_t i = start; i < uRowN; ++i, ++k)
                permuteMatrix(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteSerialDriver(std::size_t n, std::size_t r, typeVector v, bool IsRep, std::size_t uRowN,
                         std::size_t segment, std::vector<int> &z, typeMatrix &permuteMatrix) {
    
    const std::size_t first = (IsRep) ? 1u : 0u;
    auto indexMat = FromCpp14::make_unique<int[]>(segment * (r - first));
    PermuteLoadIndex(n, r, v, IsRep, segment, z, permuteMatrix, indexMat.get());
    
    std::size_t unrollRem = segment % unrollSize;
    std::size_t ind = 1u;
    std::size_t start = segment;
    std::size_t last = start + segment;
    
    for (; last <= uRowN; start += segment, last += segment, ++ind) {
        if (!IsRep) {
            typeVector vTemp(1u);
            vTemp[0u] = v[0u];
            v[0u] = v[ind];
            v[ind] = vTemp[0u];
        }
        
        PermuteWorker(r, IsRep, v, indexMat.get(), permuteMatrix, start, last, ind, first, unrollRem);
    }
    
    if (ind < static_cast<std::size_t>(v.size()) && start < uRowN) {
        const std::size_t skip = last - uRowN;
        unrollRem = uRowN % unrollSize;
        const std::size_t lastUnroll = uRowN - unrollRem;

        if (IsRep) {
            RepUnroller(r, v, indexMat.get(), permuteMatrix, start, lastUnroll, uRowN, ind);
        } else {
            typeVector vTemp(1u);
            vTemp[0u] = v[0u];
            v[0u] = v[ind];
            v[ind] = vTemp[0u];
        }

        for (std::size_t j = first, k = 0; j < r; ++j, k += skip)
            for (std::size_t i = start; i < uRowN; ++i, ++k)
                permuteMatrix(i, j) = v[indexMat[k]];
    }
}

template <typename typeMatrix, typename typeVector>
void PermuteGeneral(std::size_t n, std::size_t r, const typeVector &v, 
                    bool IsRep, std::size_t uRowN, std::vector<int> &z,
                    int intCount, typeMatrix &permuteMatrix) {
    
    const std::size_t maxInd = n - 1u;
    const std::size_t lastCol = r - 1u;
    
    if (IsRep) {
        const int maxIndInt = maxInd;
        
        for (std::size_t count = intCount; count < uRowN; ++count) {
            for (std::size_t j = 0u; j < r; ++j)
                permuteMatrix(count, j) = v[z[j]];
            
            for (int k = lastCol; k >= 0; --k) {
                if (z[k] != maxIndInt) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }
    } else {
        const std::size_t numR1 = uRowN - 1u;
        auto arrPerm = FromCpp14::make_unique<int[]>(n);
        
        for (std::size_t i = 0u; i < n; ++i)
            arrPerm[i] = z[i];
        
        if (r == n) {
            for (std::size_t count = intCount; count < numR1; ++count) {
                for (std::size_t j = 0u; j < r; ++j)
                    permuteMatrix(count, j) = v[arrPerm[j]];
                
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (std::size_t count = intCount; count < numR1; ++count) {
                for (std::size_t j = 0u; j < r; ++j)
                    permuteMatrix(count, j) = v[arrPerm[j]];
                
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0u; j < r; ++j)
            permuteMatrix(numR1, j) = v[arrPerm[j]];
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermutation(std::size_t n, std::size_t r, const typeVector &v, std::size_t uRowN,
                         std::vector<int> &z, int intCount, typeMatrix &permuteMatrix) {
    
    const std::size_t lenFreqs = z.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);
    
    const std::size_t numR1 = uRowN - 1;
    const std::size_t lastCol = r - 1;
    const std::size_t maxInd = lenFreqs - 1;
    
    for (std::size_t j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (r == lenFreqs) {
        for (std::size_t count = intCount; count < numR1; ++count) {
            for (std::size_t j = 0; j < r; ++j)
                permuteMatrix(count, j) = v[arrPerm[j]];
            
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = intCount; count < numR1; ++count) {
            for (std::size_t j = 0; j < r; ++j)
                permuteMatrix(count, j) = v[arrPerm[j]];
            
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }
    
    // Get last permutation
    for (std::size_t j = 0; j < r; ++j)
        permuteMatrix(numR1, j) = v[arrPerm[j]];
}

template <typename typeVector>
void PermutationApplyFun(std::size_t n, std::size_t r, const typeVector &v, bool IsRep,
                         std::size_t uRowN, bool Multi, std::vector<int> &z,
                         int intCount, SEXP sexpFun, SEXP rho, Rcpp::List &myList) {
    
    const std::size_t lenFreqs = (Multi) ? z.size() : 0;
    typeVector vectorPass(r);
    
    const std::size_t numR1 = uRowN - 1;
    const std::size_t lastCol = r - 1;
    const std::size_t maxInd = (Multi) ? (lenFreqs - 1) : (n - 1);
    
    if (IsRep) {
        const int maxIndInt = maxInd;
        
        for (std::size_t count = intCount; count < uRowN; ++count) {
            for (std::size_t j = 0; j < r; ++j)
                vectorPass[j] = v[z[j]];

            SETCADR(sexpFun, vectorPass);
            myList[count] = Rf_eval(sexpFun, rho);

            for (int k = lastCol; k >= 0; --k) {
                if (z[k] != maxIndInt) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }
    } else {
        const std::size_t arrLength = maxInd + 1;
        auto arrPerm = FromCpp14::make_unique<int[]>(arrLength);
        
        for (std::size_t i = 0; i < arrLength; ++i)
            arrPerm[i] = z[i];
        
        if (r == n || r == lenFreqs) {
            for (std::size_t count = intCount; count < numR1; ++count) {
                for (std::size_t j = 0; j < r; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (std::size_t count = intCount; count < numR1; ++count) {
                for (std::size_t j = 0; j < r; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                    
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < r; ++j)
            vectorPass[j] = v[arrPerm[j]];

        SETCADR(sexpFun, vectorPass);
        myList[numR1] = Rf_eval(sexpFun, rho);
    }
}

#endif
