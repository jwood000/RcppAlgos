#ifndef COMBO_RESULTS_H
#define COMBO_RESULTS_H

#include "UserConstraintFuns.h"

template <typename typeMatrix, typename typeVector>
void ComboGenResNoRep(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                      std::vector<int> z, int n, int m, int strt, int nRows,
                      const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<typeVector> vPass(m);
    
    for (int count = strt, m2 = m - 2,
         m1 = m - 1, nMinusM = n - m; count < nRows; ) {
        int numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
        }

        for (int i = m2; i >= 0; i--) {
            if (z[i] != (nMinusM + i)) {
                ++z[i];
                
                for (int j = i; j < m1; ++j)
                    z[j + 1] = z[j] + 1;

                break;
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void ComboGenResRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, 
                    std::vector<int> z, int n, int m, int strt, int nRows,
                    const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<typeVector> vPass(m);
    
    for (int count = strt, m2 = m - 2,
         m1 = m - 1, lastElement = n - 1; count < nRows; ) {
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
        }
        
        for (int i = m2; i >= 0; i--) {
            if (z[i] != lastElement) {
                ++z[i];
                
                for (int j = i; j < m1; ++j)
                    z[j + 1] = z[j];
                
                break;
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetComboResult(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                         std::vector<int> z, int n, int m, int strt, int nRows,
                         const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<int> zIndex(n), zGroup(m);
    std::vector<typeVector> vPass(m);
    
    for (int i = 0; i < n; ++i)
        zIndex[i] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    int pentExtreme = freqs.size() - m;
    
    for (int count = strt, m1 = m - 1, m2 = m - 2; count < nRows;) {
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = v[freqs[zIndex[z[j]]]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
        }
        
        for (int i = m2; i >= 0; --i) {
            if (freqs[zIndex[z[i]]] != freqs[pentExtreme + i]) {
                ++z[i];
                zGroup[i] = zIndex[z[i]];
                
                for (int j = (i + 1); j < m; ++j) {
                    zGroup[j] = zGroup[j - 1] + 1;
                    z[j] = freqs[zGroup[j]];
                }
                
                break;
            }
        }
    }
}

#endif
