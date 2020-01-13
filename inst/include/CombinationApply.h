#ifndef COMBINATION_APPLY_H
#define COMBINATION_APPLY_H

#include "NextStandard.h"

template <typename typeVector>
void ComboGeneralApplyFun(Rcpp::List &myList, const typeVector &v, std::vector<int> z,
                          int n, int m, bool IsRep, int nRows, SEXP sexpFun, SEXP rho) {
    
    typeVector vectorPass(m);
    
    if (IsRep) {
        for (int count = 0, m1 = m - 1, n1 = n - 1; count < nRows; ) {
            
            int numIter = n - z[m1];
            
            if (numIter + count > nRows)
                numIter = nRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
            }
            
            nextCombSecRep(z, m1, n1);
        }
    } else {
        for (int count = 0, m1 = m - 1, nMinusM = n - m; count < nRows;) {
            
            int numIter = n - z[m1];
            
            if ((numIter + count) > nRows)
                numIter = nRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[m1]){
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
            }
            
            nextCombSec(z, m1, nMinusM);
        }
    }
}

template <typename typeVector>
void MultisetComboApplyFun(Rcpp::List &myList, const typeVector &v, std::vector<int> z,
                           int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                           const std::vector<int> &freqs) {
    
    std::vector<int> zIndex(n);
    
    for (int i = 0; i < n; ++i)
        zIndex[i] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
    
    typeVector vectorPass(m);
    
    for (int count = 0, m1 = m - 1, 
         pentExtreme = freqs.size() - m; count < nRows;) {
        
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j)
                vectorPass[j] = v[z[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[count] = Rf_eval(sexpFun, rho);
        }
        
        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

#endif
