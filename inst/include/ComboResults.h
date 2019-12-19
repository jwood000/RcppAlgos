#ifndef COMBO_RESULTS_H
#define COMBO_RESULTS_H

#include "UserConstraintFuns.h"
#include "NextCombinatorics.h"

template <typename typeMatrix, typename typeVector>
void ComboGenResNoRep(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                      std::vector<int> z, int n, int m, int strt, int nRows,
                      const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<typeVector> vPass(m);
    
    for (int count = strt, m1 = m - 1, nMinusM = n - m; count < nRows; ) {
        
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

        nextComb(z, m1, nMinusM);
    }
}

template <typename typeMatrix, typename typeVector>
void ComboGenResRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, 
                    std::vector<int> z, int n, int m, int strt, int nRows,
                    const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<typeVector> vPass(m);
    
    for (int count = strt, m1 = m - 1, n1 = n - 1; count < nRows; ) {
        
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
        
        nextCombRep(z, m1, n1);
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
    
    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (int count = strt, m1 = m - 1, 
         pentExtreme = freqs.size() - m; count < nRows;) {
        
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
        
        nextCombMulti(freqs, zIndex, zGroup, z, m, pentExtreme);
    }
}

#endif
