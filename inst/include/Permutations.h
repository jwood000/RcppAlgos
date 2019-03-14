#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "CombPermUtils.h"
#include <memory>

template <typename typeMatrix, typename typeVector>
void PermuteGeneral(int n, int r, typeVector &v, bool repetition, int numRows,
                    std::vector<int> &z, int intCount, bool nonTrivial,
                    typeMatrix &permuteMatrix) {
    
    const unsigned long int uN = n;
    const unsigned long int uR = r;
    const unsigned long int uRowN = numRows;
    const unsigned long int lastElem = n - 1;
    const unsigned long int lastCol = r - 1;
    const unsigned long int pentultimate = n - 2;
    unsigned long int count = intCount;
    
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
        
        const unsigned long int numR1 = numRows - 1;
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
            unsigned long int phaseOne, maxN = NumPermsNoRep(n, r);
            unsigned long int segment = maxN / uN;
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
            
            unsigned long int start = segment, last = 2 * segment;
            unsigned long int ind = 1;
            typeVector vTemp(1);
            
            for (; start < uRowN; start += segment, last += segment, ++ind) {
                vTemp[0] = v[0];
                v[0] = v[ind];
                v[ind] = vTemp[0];
    
                if (last > uRowN) {
                    unsigned long int skip = last - uRowN;
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
    
    const unsigned long int lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);
    
    const unsigned long int uR = r;
    const unsigned long int uN = n;
    const unsigned long int numR1 = numRows - 1;
    const unsigned long int lastCol = r - 1;
    const unsigned long int lastElem = lenFreqs - 1;
    
    unsigned long int count = intCount;
    
    for (std::size_t j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (uR == lenFreqs) {
        const unsigned long int pentultimate = lenFreqs - 2;
        
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
    
    const unsigned long int uR = r;
    const unsigned long int uN = n;
    const unsigned long int lenFreqs = (Multi) ? z.size() : 0;
    typeVector vectorPass(r);
    
    const unsigned long int numR1 = numRows - 1;
    const unsigned long int uRowN = numRows, lastCol = uR - 1;
    const unsigned long int lastElem = (Multi) ? (lenFreqs - 1) : (n - 1);
    unsigned long int count = intCount;
    
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
        const unsigned long int arrLength = lastElem + 1;
        auto arrPerm = std::make_unique<int[]>(arrLength);
        
        for (std::size_t i = 0; i < arrLength; ++i)
            arrPerm[i] = z[i];
        
        if (uR == uN || uR == lenFreqs) {
            const unsigned long int pentultimate = lastElem - 1;
            
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
