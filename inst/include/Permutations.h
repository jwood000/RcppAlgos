#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "CombPermUtils.h"
#include <memory>

template <typename typeMatrix, typename typeVector>
void PermuteGeneral(int n, int r, typeVector &v, bool repetition, int numRows,
                    std::vector<int> &z, int intCount, bool nonTrivial,
                    typeMatrix &permuteMatrix) {
    
    const std::size_t uN = n;
    const std::size_t uR = r;
    const std::size_t uRowN = numRows;
    const std::size_t lastElem = n - 1;
    const std::size_t lastCol = r - 1;
    const std::size_t pentultimate = n - 2;
    std::size_t count = intCount;
    
    if (repetition) {
        const int lastElemInt = lastElem;
        
        for (; count < uRowN; ++count) {
            for (std::size_t j = 0; j < uR; ++j)
                permuteMatrix(count, j) = v[z[j]];
            
            for (int k = lastCol; k >= 0; --k) {
                if (z[k] != lastElemInt) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }
        
    } else if (nonTrivial) {
        
        const std::size_t numR1 = numRows - 1;
        auto arrPerm = std::make_unique<int[]>(uN);

        for (std::size_t i = 0; i < uN; ++i)
            arrPerm[i] = z[i];
        
        if (r == n) {
            for (; count < numR1; ++count) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(count, j) = v[arrPerm[j]];
                
                nextFullPerm(arrPerm.get(), lastElem, pentultimate);
            }
        } else {
            for (; count < numR1; ++count) {
                for (std::size_t j = 0; j < uR; ++j)
                    permuteMatrix(count, j) = v[arrPerm[j]];
                
                nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < uR; ++j)
            permuteMatrix(numR1, j) = v[arrPerm[j]];
        
    } else {
        
        if (n > 1) {
            std::size_t phaseOne, maxN = NumPermsNoRep(n, r);
            std::size_t segment = maxN / uN;
            phaseOne = (uRowN < segment) ? uRowN : segment;
            
            auto indexMat = std::make_unique<int[]>(phaseOne * uR);
            auto arrPerm = std::make_unique<int[]>(uN);
            
            for (int i = 0; i < n; ++i)
                arrPerm[i] = i;

            if (r == n) {
                for (std::size_t i = 0; i < phaseOne; ++i) {
                    for (std::size_t j = 0, k = i; j < uR; ++j, k += phaseOne) {
                        permuteMatrix(i, j) = v[arrPerm[j]];
                        indexMat[k] = arrPerm[j];
                    }
                    
                    nextFullPerm(arrPerm.get(), lastElem, pentultimate);
                }
            } else {
                for (std::size_t i = 0; i < phaseOne; ++i) {
                    for (std::size_t j = 0, k = i; j < uR; ++j, k += phaseOne) {
                        permuteMatrix(i, j) = v[arrPerm[j]];
                        indexMat[k] = arrPerm[j];
                    }
                    
                    nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
                }
            }
            
            std::size_t start = segment, last = 2 * segment;
            std::size_t ind = 1;
            typeVector vTemp(1);
            
            for (; start < uRowN; start += segment, last += segment, ++ind) {
                vTemp[0] = v[0];
                v[0] = v[ind];
                v[ind] = vTemp[0];
    
                if (last > uRowN) {
                    std::size_t skip = last - uRowN;
                    last = uRowN;
                    
                    for (std::size_t j = 0, k = 0; j < uR; ++j, k += skip)
                        for (std::size_t i = start; i < last; ++i, ++k)
                            permuteMatrix(i, j) = v[indexMat[k]];
                } else {
                    for (std::size_t j = 0, k = 0; j < uR; ++j)
                        for (std::size_t i = start; i < last; ++i, ++k)
                            permuteMatrix(i, j) = v[indexMat[k]];
                }
            }
        } else {
            permuteMatrix(0, 0) = v[0];
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermutation(int n, int r, typeVector &v, int numRows, std::vector<int> &z,
                         int intCount, typeMatrix &permuteMatrix) {
    
    const std::size_t lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);
    
    const std::size_t uR = r;
    const std::size_t uN = n;
    const std::size_t numR1 = numRows - 1;
    const std::size_t lastCol = r - 1;
    const std::size_t lastElem = lenFreqs - 1;
    
    std::size_t count = intCount;
    
    for (std::size_t j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (uR == lenFreqs) {
        const std::size_t pentultimate = lenFreqs - 2;
        
        for (; count < numR1; ++count) {
            for (std::size_t j = 0; j < uR; ++j)
                permuteMatrix(count, j) = v[arrPerm[j]];
            
            nextFullPerm(arrPerm.get(), lastElem, pentultimate);
        }
    } else {
        for (; count < numR1; ++count) {
            for (std::size_t j = 0; j < uR; ++j)
                permuteMatrix(count, j) = v[arrPerm[j]];
            
            nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
        }
    }
    
    // Get last permutation
    for (std::size_t j = 0; j < uR; ++j)
        permuteMatrix(numR1, j) = v[arrPerm[j]];
}

template <typename typeVector>
void PermutationApplyFun(int n, int r, typeVector &v, bool repetition,
                         int numRows, bool Multi, std::vector<int> &z,
                         int intCount, SEXP sexpFun, SEXP rho, SEXP &ans) {
    
    const std::size_t uR = r;
    const std::size_t uN = n;
    const std::size_t lenFreqs = (Multi) ? z.size() : 0;
    typeVector vectorPass(r);
    
    const std::size_t numR1 = numRows - 1;
    const std::size_t uRowN = numRows, lastCol = uR - 1;
    const std::size_t lastElem = (Multi) ? (lenFreqs - 1) : (n - 1);
    std::size_t count = intCount;
    
    if (repetition) {
        const int lastElemInt = lastElem;
        
        for (; count < uRowN; ++count) {
            for (std::size_t j = 0; j < uR; ++j)
                vectorPass[j] = v[z[j]];

            SETCADR(sexpFun, vectorPass);
            SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));

            for (int k = lastCol; k >= 0; --k) {
                if (z[k] != lastElemInt) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }
    } else {
        const std::size_t arrLength = lastElem + 1;
        auto arrPerm = std::make_unique<int[]>(arrLength);
        
        for (std::size_t i = 0; i < arrLength; ++i)
            arrPerm[i] = z[i];
        
        if (uR == uN || uR == lenFreqs) {
            const std::size_t pentultimate = lastElem - 1;
            
            for (; count < numR1; ++count) {
                for (std::size_t j = 0; j < uR; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));
                nextFullPerm(arrPerm.get(), lastElem, pentultimate);
            }
        } else {
            for (; count < numR1; ++count) {
                for (std::size_t j = 0; j < uR; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                    
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));
                nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < uR; ++j)
            vectorPass[j] = v[arrPerm[j]];

        SETCADR(sexpFun, vectorPass);
        SET_VECTOR_ELT(ans, numR1, Rf_eval(sexpFun, rho));
    }
}

#endif
