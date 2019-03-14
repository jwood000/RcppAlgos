#ifndef PERMUTE_RESULTS_H
#define PERMUTE_RESULTS_H

#include "CombPermUtils.h"
#include "ConstraintsUtils.h"
#include <memory>

template <typename typeMatrix, typename typeVector>
void PermuteGenRes(int n, int r, std::vector<typeVector> &v, bool repetition,
                   int numRows, std::vector<int> &z, unsigned long int count,
                   bool nonTrivial, typeMatrix &permuteMatrix, funcPtr<typeVector> myFun) {
    
    const unsigned long int uN = n;
    const unsigned long int uR = r;
    const unsigned long int uRowN = numRows;
    const unsigned long int lastElem = n - 1;
    const unsigned long int lastCol = r - 1;
    const unsigned long int pentultimate = n - 2;
    std::vector<typeVector> vPass(uR);
    
    if (repetition) {
        const int lastElemInt = lastElem;
        
        for (; count < uRowN; ++count) {
            for (std::size_t j = 0; j < uR; ++j) {
                vPass[j] = v[z[j]];
                permuteMatrix(count, j) = vPass[j];
            }
            
            permuteMatrix(count, uR) = myFun(vPass, uR);
            
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
                for (std::size_t j = 0; j < uR; ++j) {
                    vPass[j] = v[arrPerm[j]];
                    permuteMatrix(count, j) = vPass[j];
                }
                
                permuteMatrix(count, uR) = myFun(vPass, uR);
                nextFullPerm(arrPerm.get(), lastElem, pentultimate);
            }
        } else {
            for (; count < numR1; ++count) {
                for (std::size_t j = 0; j < uR; ++j) {
                    vPass[j] = v[arrPerm[j]];
                    permuteMatrix(count, j) = vPass[j];
                }
                
                permuteMatrix(count, uR) = myFun(vPass, uR);
                nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0; j < uR; ++j) {
            vPass[j] = v[arrPerm[j]];
            permuteMatrix(numR1, j) = vPass[j];
        }
        
        permuteMatrix(numR1, uR) = myFun(vPass, uR);
        
    } else {
        
        if (n > 1) {
            unsigned long int phaseOne, maxN = NumPermsNoRep(n, r);
            unsigned long int segment = maxN / uN;
            phaseOne = (uRowN < segment) ? uRowN : segment;

            auto indexMat = std::make_unique<int[]>(phaseOne * uR);
            auto arrPerm = std::make_unique<int[]>(uN);
            
            for (std::size_t i = 0; i < uN; ++i)
                arrPerm[i] = static_cast<int>(i);

            if (r == n) {
                for (std::size_t i = 0, k = 0; i < phaseOne; ++i) {
                    for (std::size_t j = 0; j < uR; ++j, ++k) {
                        vPass[j] = v[arrPerm[j]];
                        permuteMatrix(i, j) = vPass[j];
                        indexMat[k] = arrPerm[j];
                    }

                    permuteMatrix(i, uR) = myFun(vPass, uR);
                    nextFullPerm(arrPerm.get(), lastElem, pentultimate);
                }
            } else {
                for (std::size_t i = 0, k = 0; i < phaseOne; ++i) {
                    for (std::size_t j = 0; j < uR; ++j, ++k) {
                        vPass[j] = v[arrPerm[j]];
                        permuteMatrix(i, j) = vPass[j];
                        indexMat[k] = arrPerm[j];
                    }
                    
                    permuteMatrix(i, uR) = myFun(vPass, uR);
                    nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
                }
            }
            
            unsigned long int start = segment, last = 2 * segment;
            unsigned long int ind = 1;
            std::vector<typeVector> vTemp(1);
            
            for (; start < uRowN; start += segment, last += segment, ++ind) {
                vTemp[0] = v[0];
                v[0] = v[ind];
                v[ind] = vTemp[0];
    
                if (last > uRowN)
                    last = uRowN;
                
                for (std::size_t i = start, k = 0; i < last; ++i) {
                    for (std::size_t j = 0; j < uR; ++j, ++k) {
                        vPass[j] = v[indexMat[k]];
                        permuteMatrix(i, j) = vPass[j];
                    }
                    
                    permuteMatrix(i, uR) = myFun(vPass, uR);
                }
            }
        } else {
            permuteMatrix(0, 0) = v[0];
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermRes(int n, int r, std::vector<typeVector> &v, int numRows,
                     unsigned long int count, std::vector<int> &z,
                     typeMatrix &permuteMatrix, funcPtr<typeVector> myFun) {
    
    const unsigned long int lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);
    std::vector<typeVector> vPass(r);
    
    const unsigned long int uR = r;
    const unsigned long int uN = n;
    const unsigned long int lastCol = r - 1;
    const unsigned long int numR1 = numRows - 1;
    const unsigned long int lastElem = lenFreqs - 1;
    
    for (std::size_t j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (uR == lenFreqs) {
        const unsigned long int pentultimate = lenFreqs - 2;
        
        for (; count < numR1; ++count) {
            for (std::size_t j = 0; j < uR; ++j) {
                vPass[j] = v[arrPerm[j]];
                permuteMatrix(count, j) = vPass[j];
            }
            
            permuteMatrix(count, uR) = myFun(vPass, uR);
            nextFullPerm(arrPerm.get(), lastElem, pentultimate);
        }
    } else {
        for (; count < numR1; ++count) {
            for (std::size_t j = 0; j < uR; ++j) {
                vPass[j] = v[arrPerm[j]];
                permuteMatrix(count, j) = vPass[j];
            }
            
            permuteMatrix(count, uR) = myFun(vPass, uR);
            nextPartialPerm(arrPerm.get(), uR, lastCol, uN, lastElem);
        }
    }
    
    // Get last permutation
    for (std::size_t j = 0; j < uR; ++j) {
        vPass[j] = v[arrPerm[j]];
        permuteMatrix(numR1, j) = vPass[j];
    }
    
    permuteMatrix(numR1, uR) = myFun(vPass, uR);
}

#endif
