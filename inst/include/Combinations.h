#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <CombPermUtils.h>

namespace Combinations {

    template <typename typeMatrix, typename typeVector>
    typeMatrix ComboGeneral(int n, int r, typeVector &v, bool repetition, 
                            int numRows, bool xtraCol, std::vector<int> &z) {
        
        int r1 = r - 1, r2 = r - 2, numIter;
        int numCols, count = 0;
        numCols = xtraCol ? (r + 1) : r;
        typeMatrix combinationMatrix = Rcpp::no_init_matrix(numRows, numCols);
        
        if (repetition) {
            int lastElement = n - 1;
            
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
            int nMinusR = n - r;
            
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

        return combinationMatrix;
    }
    
    template <typename typeMatrix, typename typeVector>
    typeMatrix MultisetCombination(int n, int r, typeVector &v, std::vector<int> &Reps,
                                   std::vector<int> &freqs, int numRows, 
                                   bool xtraCol, std::vector<int> &z) {
        
        std::vector<int> zIndex(n), zGroup(r);
        int numIter, count = 0, sizeFreqs = 0;
        int r1 = r - 1, r2 = r - 2;
        
        for (int i = 0; i < n; ++i) {
            zIndex[i] = sizeFreqs;
            sizeFreqs += Reps[i];
        }
        
        int numCols = xtraCol ? (r + 1) : r;
        typeMatrix combinationMatrix = Rcpp::no_init_matrix(numRows, numCols);
        
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
        
        return combinationMatrix;
    }
    
    template <typename typeVector>
    SEXP ComboGeneralApplyFun(int n, int r, typeVector &v, bool repetition, 
                              int numRows, std::vector<int> &z, SEXP func, SEXP rho) {
        
        int r1 = r - 1, r2 = r - 2;
        int numIter, count = 0;;
        
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, numRows));
        SEXP R_fcall = PROTECT(Rf_lang2(func, R_NilValue));
        typeVector vectorPass(r);
        
        if (repetition) {
            int lastElement = n - 1;
            
            while (count < numRows) {
                numIter = n - z[r1];
                
                if (numIter + count > numRows)
                    numIter = numRows - count;
                
                for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
                    for (int k = 0; k < r; ++k)
                        vectorPass[k] = v[z[k]];
                    
                    SETCADR(R_fcall, vectorPass);
                    SET_VECTOR_ELT(ans, count, Rf_eval(R_fcall, rho));
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
            
            while (count < numRows) {
                numIter = n - z[r1];
                
                if ((numIter + count) > numRows)
                    numIter = numRows - count;
                
                for (int i = 0; i < numIter; ++i, ++count, ++z[r1]){
                    for (int k = 0; k < r; ++k)
                        vectorPass[k] = v[z[k]];
                
                    SETCADR(R_fcall, vectorPass);
                    SET_VECTOR_ELT(ans, count, Rf_eval(R_fcall, rho));
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
        
        UNPROTECT(2);
        return ans;
    }
    
    template <typename typeVector>
    SEXP MultisetComboApplyFun(int n, int r, typeVector &v, std::vector<int> Reps,
                               std::vector<int> freqs, int numRows, std::vector<int> &z,
                               SEXP func, SEXP rho) {

        int count = 0, sizeFreqs = 0, numIter;

        std::vector<int> zIndex(n), zGroup(r);
        int r1 = r - 1, r2 = r - 2;
        
        for (int i = 0; i < n; ++i) {
            zIndex[i] = sizeFreqs;
            sizeFreqs += Reps[i];
        }
        
        // location in freqs that represents the maximal
        // value of the second to the last element
        int pentExtreme = sizeFreqs - r;
        
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, numRows));
        SEXP R_fcall = PROTECT(Rf_lang2(func, R_NilValue));
        typeVector vectorPass(r);

        while (count < numRows) {
            numIter = n - z[r1];

            if (numIter + count > numRows)
                numIter = numRows - count;

            for (int i = 0; i < numIter; ++i, ++count, ++z[r1]) {
                for (int k = 0; k < r; ++k)
                    vectorPass[k] = v[freqs[zIndex[z[k]]]];
            
                SETCADR(R_fcall, vectorPass);
                SET_VECTOR_ELT(ans, count, Rf_eval(R_fcall, rho));
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
        
        UNPROTECT(2);
        return ans;
    }
}

#endif
