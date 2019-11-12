#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <Rcpp.h>

template <typename typeMatrix, typename typeVector>
void CombinationsNoRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                       int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    for (int count = strt, m2 = m - 2,
         m1 = m - 1, nMinusM = n - m; count < nRows;) {
        
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1])
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[z[j]];
        
        for (int i = m2; i >= 0; --i) {
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
void CombinationsRep(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                     int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    for (int count = strt, m2 = m - 2,
         m1 = m - 1, lastElement = n - 1; count < nRows;) {
        
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1])
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[z[j]];
        
        for (int i = m2; i >= 0; --i) {
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
void MultisetCombination(typeMatrix &matRcpp, const typeVector &v, std::vector<int> z,
                         int n, int m, int strt, int nRows, const std::vector<int> &freqs) {
    
    std::vector<int> zIndex(n), zGroup(m);
    
    for (int i = 0; i < n; ++i)
        zIndex[i] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    int pentExtreme = freqs.size() - m;
    
    for (int count = strt, m1 = m - 1, m2 = m - 2; count < nRows;) {
        
        int numIter = n - z[m1];
        
        if (numIter + count > nRows)
            numIter = nRows - count;
        
        for (int i = 0; i < numIter; ++i, ++count, ++z[m1])
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[freqs[zIndex[z[j]]]];
        
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

template <typename typeVector>
void ComboGeneralApplyFun(Rcpp::List &myList, const typeVector &v, std::vector<int> &z,
                          int n, int m, bool IsRep, int nRows, SEXP sexpFun, SEXP rho) {
    
    typeVector vectorPass(m);
    
    if (IsRep) {
        for (int count = 0, m1 = m - 1, 
             m2 = m - 2, lastElement = n - 1; count < nRows; ) {
            
            int numIter = n - z[m1];
            
            if (numIter + count > nRows)
                numIter = nRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
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
    } else {
        for (int count = 0, m1 = m - 1, 
             m2 = m - 2, nMinusM = n - m; count < nRows;) {
            
            int numIter = n - z[m1];
            
            if ((numIter + count) > nRows)
                numIter = nRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count, ++z[m1]){
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
            
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
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
}

template <typename typeVector>
void MultisetComboApplyFun(Rcpp::List &myList, const typeVector &v, std::vector<int> &z,
                           int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                           const std::vector<int> &freqs) {

    std::vector<int> zIndex(n), zGroup(m);
    
    for (int i = 0; i < n; ++i)
        zIndex[i] = std::find(freqs.cbegin(), freqs.cend(), i) - freqs.cbegin();
    
    // location in freqs that represents the maximal
    // value of the second to the last element
    const int pentExtreme = freqs.size() - m;
    typeVector vectorPass(m);

    for (int count = 0, m1 = m - 1, m2 = m - 2; count < nRows;) {
        
        int numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (int i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j)
                vectorPass[j] = v[freqs[zIndex[z[j]]]];
        
            SETCADR(sexpFun, vectorPass);
            myList[count] = Rf_eval(sexpFun, rho);
        }

        for (int i = m2; i >= 0; i--) {
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
