#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <Rcpp.h>

template <typename typeMatrix, typename typeVector>
void ComboGeneral(int n, int r, typeVector &v, bool repetition, int count,
                  int numRows, std::vector<int> &z, typeMatrix &combinationMatrix) {
    
    const int r1 = r - 1;
    const int r2 = r - 2;
    int numIter;
    
    if (repetition) {
        const int lastElement = n - 1;
        
        while (count < numRows) {
            numIter = n - z[r1];
            
            if (numIter + count > numRows)
                numIter = numRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[r1])
                for (int k = 0; k < r; ++k)
                    combinationMatrix(count, k) = v[z[k]];
            
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
        const int nMinusR = n - r;
        
        while (count < numRows) {
            numIter = n - z[r1];
            
            if (numIter + count > numRows)
                numIter = numRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[r1])
                for (int k = 0; k < r; ++k)
                    combinationMatrix(count, k) = v[z[k]];
            
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
void MultisetCombination(int n, int r, typeVector &v, std::vector<int> &Reps,
                               std::vector<int> &freqs, int count, int numRows,
                               std::vector<int> &z, typeMatrix &combinationMatrix) {
    
    std::vector<int> zIndex(n), zGroup(r);
    int numIter, sizeFreqs = 0;
    const int r1 = r - 1;
    const int r2 = r - 2;
    
    for (int i = 0; i < n; ++i) {
        zIndex[i] = sizeFreqs;
        sizeFreqs += Reps[i];
    }
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    int pentExtreme = sizeFreqs - r;
    
    while (count < numRows) {
        numIter = n - z[r1];
        
        if (numIter + count > numRows)
            numIter = numRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[r1])
            for (int k = 0; k < r; ++k)
                combinationMatrix(count, k) = v[freqs[zIndex[z[k]]]];
        
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

template <typename typeVector>
void ComboGeneralApplyFun(int n, int r, typeVector &v, bool repetition, int count,
                          int numRows, std::vector<int> &z, SEXP sexpFun, SEXP rho, SEXP &ans) {
    
    const int r1 = r - 1;
    const int r2 = r - 2;
    int numIter;
    typeVector vectorPass(r);
    
    if (repetition) {
        const int lastElement = n - 1;
        
        while (count < numRows) {
            numIter = n - z[r1];
            
            if (numIter + count > numRows)
                numIter = numRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
                for (int k = 0; k < r; ++k)
                    vectorPass[k] = v[z[k]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));
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
        const int nMinusR = n - r;
        
        while (count < numRows) {
            numIter = n - z[r1];
            
            if ((numIter + count) > numRows)
                numIter = numRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[r1]){
                for (int k = 0; k < r; ++k)
                    vectorPass[k] = v[z[k]];
            
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));
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

template <typename typeVector>
void MultisetComboApplyFun(int n, int r, typeVector &v, std::vector<int> &Reps,
                           std::vector<int> &freqs, int numRows, std::vector<int> &z,
                           int count, SEXP sexpFun, SEXP rho, SEXP &ans) {

    int sizeFreqs = 0, numIter;
    std::vector<int> zIndex(n), zGroup(r);
    const int r1 = r - 1;
    const int r2 = r - 2;
    
    for (int i = 0; i < n; ++i) {
        zIndex[i] = sizeFreqs;
        sizeFreqs += Reps[i];
    }
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    const int pentExtreme = sizeFreqs - r;
    typeVector vectorPass(r);

    while (count < numRows) {
        numIter = n - z[r1];

        if (numIter + count > numRows)
            numIter = numRows - count;

        for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
            for (int k = 0; k < r; ++k)
                vectorPass[k] = v[freqs[zIndex[z[k]]]];
        
            SETCADR(sexpFun, vectorPass);
            SET_VECTOR_ELT(ans, count, Rf_eval(sexpFun, rho));
        }

        for (int i = r2; i >= 0; i--) {
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
