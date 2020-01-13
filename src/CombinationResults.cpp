#include "UserConstraintFuns.h"
#include "NextStandard.h"
#include "RMatrix.h"

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
        
        nextCombSec(z, m1, nMinusM);
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
        
        nextCombSecRep(z, m1, n1);
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetComboResult(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                         std::vector<int> z, int n, int m, int strt, int nRows,
                         const std::vector<int> &freqs, funcPtr<typeVector> myFun) {
    
    std::vector<int> zIndex(n);
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
                vPass[j] = v[z[j]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
        }
        
        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

template void ComboGenResNoRep(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                               int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void ComboGenResNoRep(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                               int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void ComboGenResNoRep(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                               int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void ComboGenResNoRep(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                               int, int, int, int, const std::vector<int>&, funcPtr<double>);

template void ComboGenResRep(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void ComboGenResRep(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void ComboGenResRep(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void ComboGenResRep(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<double>);

template void MultisetComboResult(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void MultisetComboResult(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void MultisetComboResult(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void MultisetComboResult(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<double>);
