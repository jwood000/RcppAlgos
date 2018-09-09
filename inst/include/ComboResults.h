#ifndef COMBO_RESULTS_H
#define COMBO_RESULTS_H

#include <Rcpp.h>
#include <ConstraintsUtils.h>

template <typename typeMatrix, typename typeVector>
void ComboGenRes(int n, int r, std::vector<typeVector> &v, bool repetition,
                  int nRows, int count, std::vector<int> &z,
                  typeMatrix combinationMatrix, funcPtr<typeVector> myFun) {
    
    int r1 = r - 1, r2 = r - 2, numIter;
    std::vector<typeVector> vPass(r);
    unsigned long int uR = r;
    
    if (repetition) {
        int lastElement = n - 1;
        
        while (count < nRows) {
            numIter = n - z[r1];
            
            if (numIter + count > nRows)
                numIter = nRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
                for (int k = 0; k < r; ++k) {
                    vPass[k] = v[z[k]];
                    combinationMatrix(count, k) = vPass[k];
                }
                
                combinationMatrix(count, r) = myFun(vPass, uR);
            }
            
            for (int i = r2; i >= 0; i--) {
                if (z[i] != lastElement) {
                    ++z[i];
                    for (int k = i; k < r1; ++k)
                        z[k + 1] = z[k];
                    
                    break;
                }
            }
        }
    } else {
        int nMinusR = n - r;
        
        while (count < nRows) {
            numIter = n - z[r1];

            if (numIter + count > nRows)
                numIter = nRows - count;

            for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
                for (int k = 0; k < r; ++k) {
                    vPass[k] = v[z[k]];
                    combinationMatrix(count, k) = vPass[k];
                }
                
                combinationMatrix(count, r) = myFun(vPass, uR);
            }

            for (int i = r2; i >= 0; i--) {
                if (z[i] != (nMinusR + i)) {
                    ++z[i];
                    for (int k = i; k < r1; ++k)
                        z[k + 1] = z[k] + 1;

                    break;
                }
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetComboResult(int n, int r, std::vector<typeVector> &v, std::vector<int> &Reps,
                         std::vector<int> &freqs, int nRows, int count,
                         std::vector<int> &z, typeMatrix combinationMatrix,
                         funcPtr<typeVector> myFun) {
    
    std::vector<int> zIndex(n), zGroup(r);
    std::vector<typeVector> vPass(r);
    int numIter, sizeFreqs = 0;
    int r1 = r - 1, r2 = r - 2;
    unsigned long int uR = r;
    
    for (int i = 0; i < n; ++i) {
        zIndex[i] = sizeFreqs;
        sizeFreqs += Reps[i];
    }
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    int pentExtreme = sizeFreqs - r;
    
    while (count < nRows) {
        numIter = n - z[r1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
            for (int k = 0; k < r; ++k) {
                vPass[k] = v[freqs[zIndex[z[k]]]];
                combinationMatrix(count, k) = vPass[k];
            }
            
            combinationMatrix(count, r) = myFun(vPass, uR);
        }
        
        for (int i = r2; i >= 0; --i) {
            if (freqs[zIndex[z[i]]] != freqs[pentExtreme + i]) {
                ++z[i];
                zGroup[i] = zIndex[z[i]];
                
                for (int k = (i + 1); k < r; ++k) {
                    zGroup[k] = zGroup[k - 1] + 1;
                    z[k] = freqs[zGroup[k]];
                }
                
                break;
            }
        }
    }
}

#endif
