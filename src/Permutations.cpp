#include <Permutations.h>
#include <ConstraintsUtils.h>
using namespace Rcpp;

SEXP PermutationsRcpp(int n, int m, bool repetition, CharacterVector vStr,
                      int nRows, std::vector<int> vInt, std::vector<double> vNum,
                      bool isMult, bool isFac, bool keepRes, std::vector<int> startZ,
                      bool isChar, SEXP Rv, bool isInt, std::vector<int> myReps, 
                      SEXP f1, SEXP f2, bool nonTrivial) {
    
    if (isChar) {
        if (isMult) {
            return Permutations::MultisetPermutation<CharacterMatrix>(n, m, vStr, myReps, nRows, false, startZ);
        } else {
            return Permutations::PermuteGeneral<CharacterMatrix>(n, m, vStr, repetition, 
                                                                 nRows, false, startZ, nonTrivial);
        }
    } else {
        if (Rf_isNull(f1))
            keepRes = false;
        
        if (keepRes) {
            NumericMatrix matRes;
            
            std::string mainFun2 = as<std::string >(f1);
            if (mainFun2 != "prod" && mainFun2 != "sum" && mainFun2 != "mean"
                    && mainFun2 != "max" && mainFun2 != "min") {
                stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            std::vector<double> rowVec(m);
            XPtr<funcPtr> xpFun2 = putFunPtrInXPtr(mainFun2);
            funcPtr myFun2 = *xpFun2;
            
            if (isMult) {
                matRes = Permutations::MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, true, startZ);
            } else {
                matRes = Permutations::PermuteGeneral<NumericMatrix>(n, m, vNum, repetition,
                                                                     nRows, true, startZ, nonTrivial);
            }
            
            for (std::size_t i = 0; i < nRows; i++) {
                for (std::size_t j = 0; j < m; j++)
                    rowVec[j] = matRes(i, j);
                
                matRes(i, m) = myFun2(rowVec);
            }
            
            return matRes;
        } else {
            if (isFac) {
                IntegerMatrix factorMat;
                IntegerVector testFactor = as<IntegerVector>(Rv);
                CharacterVector myClass = testFactor.attr("class");
                CharacterVector myLevels = testFactor.attr("levels");
                
                if (isMult) {
                    factorMat = Permutations::MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false, startZ);
                } else {
                    factorMat = Permutations::PermuteGeneral<IntegerMatrix>(n, m, vInt, repetition, 
                                                                            nRows, false, startZ, nonTrivial);
                }
                
                factorMat.attr("class") = myClass;
                factorMat.attr("levels") = myLevels;
                
                return factorMat;
            } else {
                if (isInt) {
                    if (isMult) {
                        return Permutations::MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false, startZ);
                    } else {
                        return Permutations::PermuteGeneral<IntegerMatrix>(n, m, vInt, repetition, 
                                                                           nRows, false, startZ, nonTrivial);
                    }
                } else {
                    if (isMult) {
                        return Permutations::MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, false, startZ);
                    } else {
                        return Permutations::PermuteGeneral<NumericMatrix>(n, m, vNum, repetition, 
                                                                           nRows, false, startZ, nonTrivial);
                    }
                }
            }
        }
    }
}
