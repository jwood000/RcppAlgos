#ifndef PERMUTATIONS_APPLY_H
#define PERMUTATIONS_APPLY_H

#include "NextStandard.h"
#include "Cpp14MakeUnique.h"

template <int RTYPE>
void PermutationApplyFun(Rcpp::List &myList, const Rcpp::Vector<RTYPE> &v, std::vector<int> z, int n,
                         int m, bool IsRep, bool IsMult, int nRows, SEXP sexpFun, SEXP rho) {
    
    const int lenFreqs = (IsMult) ? z.size() : 0;
    Rcpp::Vector<RTYPE> vectorPass(m);
    
    const int numR1 = nRows - 1;
    const int lastCol = m - 1;
    const int maxInd = (IsMult) ? (lenFreqs - 1) : (n - 1);
    
    if (IsRep) {
        for (int count = 0; count < nRows; ++count) {
            for (int j = 0; j < m; ++j)
                vectorPass[j] = v[z[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[count] = Rf_eval(sexpFun, rho);
            
            for (int i = lastCol; i >= 0; --i) {
                if (z[i] != maxInd) {
                    ++z[i];
                    break;
                } else {
                    z[i] = 0;
                }
            }
        }
    } else {
        const int arrLength = maxInd + 1;
        auto arrPerm = FromCpp14::make_unique<int[]>(arrLength);
        
        for (int i = 0; i < arrLength; ++i)
            arrPerm[i] = z[i];
        
        if (m == n || m == lenFreqs) {
            for (int count = 0; count < numR1; ++count) {
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (int count = 0; count < numR1; ++count) {
                for (int j = 0; j < m; ++j)
                    vectorPass[j] = v[arrPerm[j]];
                
                SETCADR(sexpFun, vectorPass);
                myList[count] = Rf_eval(sexpFun, rho);
                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
        
        // Get last permutation
        for (int j = 0; j < m; ++j)
            vectorPass[j] = v[arrPerm[j]];
        
        SETCADR(sexpFun, vectorPass);
        myList[numR1] = Rf_eval(sexpFun, rho);
    }
}

#endif
