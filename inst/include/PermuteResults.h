#ifndef PERMUTE_RESULTS_H
#define PERMUTE_RESULTS_H

#include "CombPermUtils.h"
#include "ConstraintsUtils.h"
#include <memory>

template <typename typeMatrix, typename typeVector>
void PermuteGenRes(std::size_t n, std::size_t r, const std::vector<typeVector> &v, 
                   bool repetition, std::size_t uRowN, std::vector<int> &z, int intCount,
                   typeMatrix &permuteMatrix, funcPtr<typeVector> myFun) {
    
    const std::size_t maxInd = n - 1u;
    const std::size_t lastCol = r - 1u;
    std::vector<typeVector> vPass(r);
    
    if (repetition) {
        const int maxIndInt = maxInd;
        
        for (std::size_t count = intCount; count < uRowN; ++count) {
            for (std::size_t j = 0u; j < r; ++j) {
                vPass[j] = v[z[j]];
                permuteMatrix(count, j) = vPass[j];
            }
            
            permuteMatrix(count, r) = myFun(vPass, r);
            
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
        auto arrPerm = std::make_unique<int[]>(n);

        for (std::size_t i = 0u; i < n; ++i)
            arrPerm[i] = z[i];
        
        if (r == n) {
            // Since we are getting all permutations of v, we know that
            // the result of myFun on v will remain the same for the 5
            // functions defined in ConstraintsUtils.h (i.e. order does
            // not matter for min, max, prod, mean, & sum).
            for (std::size_t j = 0u; j < r; ++j) {
                vPass[j] = v[arrPerm[j]];
                permuteMatrix(intCount, j) = vPass[j];
            }
            
            const auto myRes = myFun(vPass, r);
            permuteMatrix(intCount, r) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
            
            for (std::size_t count = intCount + 1; count < numR1; ++count) {
                for (std::size_t j = 0u; j < r; ++j)
                    permuteMatrix(count, j) = v[arrPerm[j]];
                
                permuteMatrix(count, r) = myRes;
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (std::size_t count = intCount; count < numR1; ++count) {
                for (std::size_t j = 0u; j < r; ++j) {
                    vPass[j] = v[arrPerm[j]];
                    permuteMatrix(count, j) = vPass[j];
                }
                
                permuteMatrix(count, r) = myFun(vPass, r);
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
        
        // Get last permutation
        for (std::size_t j = 0u; j < r; ++j) {
            vPass[j] = v[arrPerm[j]];
            permuteMatrix(numR1, j) = vPass[j];
        }
        
        permuteMatrix(numR1, r) = myFun(vPass, r);
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermRes(std::size_t n, std::size_t r, const std::vector<typeVector> &v, 
                     std::size_t numRows, int intCount, std::vector<int> &z,
                     typeMatrix &permuteMatrix, funcPtr<typeVector> myFun) {
    
    const std::size_t lenFreqs = z.size();
    const std::size_t lastCol = r - 1u;
    auto arrPerm = std::make_unique<int[]>(lenFreqs);
    std::vector<typeVector> vPass(r);
    
    const std::size_t numR1 = numRows - 1u;
    const std::size_t maxInd = lenFreqs - 1u;
    
    for (std::size_t j = 0u; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (r == lenFreqs) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (std::size_t j = 0u; j < r; ++j) {
            vPass[j] = v[arrPerm[j]];
            permuteMatrix(intCount, j) = vPass[j];
        }
        
        const auto myRes = myFun(vPass, r);
        permuteMatrix(intCount, r) = myRes;
        nextFullPerm(arrPerm.get(), maxInd);
        
        for (std::size_t count = intCount + 1; count < numR1; ++count) {
            for (std::size_t j = 0u; j < r; ++j)
                permuteMatrix(count, j) = v[arrPerm[j]];
            
            permuteMatrix(count, r) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = intCount; count < numR1; ++count) {
            for (std::size_t j = 0u; j < r; ++j) {
                vPass[j] = v[arrPerm[j]];
                permuteMatrix(count, j) = vPass[j];
            }
            
            permuteMatrix(count, r) = myFun(vPass, r);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }
    
    // Get last permutation
    for (std::size_t j = 0u; j < r; ++j) {
        vPass[j] = v[arrPerm[j]];
        permuteMatrix(numR1, j) = vPass[j];
    }
    
    permuteMatrix(numR1, r) = myFun(vPass, r);
}

#endif
