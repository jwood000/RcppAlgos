#include "NextStandard.h"
#include "Cpp14MakeUnique.h"
#include "UserConstraintFuns.h"
#include "RMatrix.h"

template <typename typeMatrix, typename typeVector>
void PermuteGenResNoRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                        int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                        funcPtr<typeVector> myFun) {
    
    const int maxInd = n - 1;
    const int numR1 = nRows - 1;
    
    std::vector<typeVector> vPass(m);
    auto arrPerm = FromCpp14::make_unique<int[]>(n);
    
    for (int i = 0; i < n; ++i)
        arrPerm[i] = z[i];
    
    if (m == n) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            matRcpp(strt, j) = vPass[j];
        }
        
        const auto myRes = myFun(vPass, m);
        matRcpp(strt, m) = myRes;
        nextFullPerm(arrPerm.get(), maxInd);
        
        for (int count = strt + 1; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];
            
            matRcpp(count, m) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        const int lastCol = m - 1;
        
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }
    
    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        matRcpp(numR1, j) = vPass[j];
    }
    
    matRcpp(numR1, m) = myFun(vPass, m);
}

template <typename typeMatrix, typename typeVector>
void PermuteGenResRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                      int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                      funcPtr<typeVector> myFun) {
    
    const int maxInd = n - 1;
    const int lastCol = m - 1;
    std::vector<typeVector> vPass(m);
    
    for (int count = strt; count < nRows; ++count) {
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[z[j]];
            matRcpp(count, j) = vPass[j];
        }
        
        matRcpp(count, m) = myFun(vPass, m);
        
        for (int i = lastCol; i >= 0; --i) {
            if (z[i] != maxInd) {
                ++z[i];
                break;
            } else {
                z[i] = 0;
            }
        }
    }
}

template <typename typeMatrix, typename typeVector>
void MultisetPermRes(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                     int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                     funcPtr<typeVector> myFun) {
    
    const int lenFreqs = freqs.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);
    std::vector<typeVector> vPass(m);
    
    const int numR1 = nRows - 1;
    const int maxInd = lenFreqs - 1;
    
    for (int j = 0; j < lenFreqs; ++j)
        arrPerm[j] = z[j];
    
    if (m == lenFreqs) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            matRcpp(strt, j) = vPass[j];
        }
        
        const auto myRes = myFun(vPass, m);
        matRcpp(strt, m) = myRes;
        nextFullPerm(arrPerm.get(), maxInd);
        
        for (int count = strt + 1; count < numR1; ++count) {
            for (int j = 0; j < m; ++j)
                matRcpp(count, j) = v[arrPerm[j]];
            
            matRcpp(count, m) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        const int lastCol = m - 1;
        
        for (int count = strt; count < numR1; ++count) {
            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                matRcpp(count, j) = vPass[j];
            }
            
            matRcpp(count, m) = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }
    
    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        matRcpp(numR1, j) = vPass[j];
    }
    
    matRcpp(numR1, m) = myFun(vPass, m);
}

template void PermuteGenResNoRep(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                                 int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void PermuteGenResNoRep(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                                 int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void PermuteGenResNoRep(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                                 int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void PermuteGenResNoRep(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                                 int, int, int, int, const std::vector<int>&, funcPtr<double>);

template void PermuteGenResRep(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void PermuteGenResRep(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void PermuteGenResRep(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void PermuteGenResRep(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                             int, int, int, int, const std::vector<int>&, funcPtr<double>);

template void MultisetPermRes(Rcpp::IntegerMatrix&, const std::vector<int>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void MultisetPermRes(Rcpp::NumericMatrix&, const std::vector<double>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<double>);
template void MultisetPermRes(RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<int>);
template void MultisetPermRes(RcppParallel::RMatrix<double>&, const std::vector<double>&, std::vector<int>,
                                  int, int, int, int, const std::vector<int>&, funcPtr<double>);

