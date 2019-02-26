#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include "Combinations.h"
#include "Permutations.h"
#include "ConstraintsUtils.h"
#include "ComboResults.h"
#include "PermuteResults.h"
#include "NthResult.h"
#include "CountGmp.h"

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

template <typename typeRcpp>
typeRcpp SubMat(typeRcpp m, int n) {
    int k = m.ncol();
    typeRcpp subMatrix = Rcpp::no_init_matrix(n, k);
    
    for (int i = 0; i < n; ++i)
        subMatrix(i, Rcpp::_) = m(i, Rcpp::_);
    
    return subMatrix;
}

void CharacterReturn(int n, int m, Rcpp::CharacterVector v, bool IsRep, int nRows,
                     bool IsComb, std::vector<int> myReps, std::vector<int> freqs,
                     std::vector<int> z, bool permNonTriv, bool IsMultiset,
                     bool keepRes, Rcpp::CharacterMatrix &matRcpp, int count) {
    if (IsComb) {
        if (IsMultiset)
            MultisetCombination(n, m, v, myReps, freqs, count, nRows, z, matRcpp);
        else
            ComboGeneral(n, m, v, IsRep, count, nRows, z, matRcpp);
    } else {
        if (IsMultiset)
            MultisetPermutation(n, m, v, nRows, z, count, matRcpp);
        else
            PermuteGeneral(n, m, v, IsRep, nRows, z, count, permNonTriv, matRcpp);
    }
}

template <typename typeRcpp, typename typeVector>
void GeneralReturn(int n, int m, std::vector<typeVector> v, bool IsRep, int nRows, bool IsComb,
                   std::vector<int> myReps, std::vector<int> freqs, std::vector<int> z,
                   bool permNonTriv, bool IsMultiset, funcPtr<typeVector> myFun,
                   bool keepRes, typeRcpp &matRcpp, int count) {
    if (keepRes) {
        if (IsComb) {
            if (IsMultiset)
                MultisetComboResult(n, m, v, myReps, freqs, nRows, count, z, matRcpp, myFun);
            else
                ComboGenRes(n, m, v, IsRep, nRows, count, z, matRcpp, myFun);
        } else {
            if (IsMultiset)
                MultisetPermRes(n, m, v, nRows, count, z, matRcpp, myFun);
            else
                PermuteGenRes(n, m, v, IsRep, nRows, z, count, permNonTriv, matRcpp, myFun);
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                MultisetCombination(n, m, v, myReps, freqs, count, nRows, z, matRcpp);
            else
                ComboGeneral(n, m, v, IsRep, count, nRows, z, matRcpp);
        } else {
            if (IsMultiset)
                MultisetPermutation(n, m, v, nRows, z, count, matRcpp);
            else
                PermuteGeneral(n, m, v, IsRep, nRows, z, count, permNonTriv, matRcpp);
        }
    }
}

// This is called when we can't easily produce a monotonic sequence overall,
// and we must generate and test every possible combination/permutation
template <typename typeRcpp, typename typeVector>
typeRcpp SpecCaseRet(int n, int m, std::vector<typeVector> v, bool IsRep, int nRows, 
                     bool keepRes, std::vector<int> z, double lower, std::string mainFun,
                     bool IsMultiset, double computedRows, std::vector<std::string> compFunVec,
                     std::vector<typeVector> myLim, bool IsComb, std::vector<int> myReps,
                     std::vector<int> freqs, bool bLower, bool permNonTriv, double userRows) {
    
    if (!bLower) {
        if (computedRows > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
        
        nRows = static_cast<int>(computedRows);
    }
    
    std::vector<typeVector> rowVec(m);
    bool Success = false;
    std::vector<int> indexMatch;
    indexMatch.reserve(nRows);
    
    Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(mainFun);
    funcPtr<typeVector> myFun = *xpFun;
    typeRcpp matRes = Rcpp::no_init_matrix(nRows, m + 1);
    typeVector testVal;
    
    GeneralReturn(n, m, v, IsRep, nRows, IsComb, myReps, freqs, 
                  z, permNonTriv, IsMultiset, myFun, true, matRes, 0);
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp = putCompPtrInXPtr<typeVector>(compFunVec[0]);
    compPtr<typeVector> myComp = *xpComp;
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp2 = xpComp;
    compPtr<typeVector> myComp2;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            testVal = matRes(i, m);
            Success = myComp(testVal, myLim);
            
            if (Success)
                indexMatch.push_back(i);
        }
    } else {
        xpComp2 = putCompPtrInXPtr<typeVector>(compFunVec[1]);
        myComp2 = *xpComp2;
        std::vector<typeVector> myLim2 = myLim;
        myLim2.erase(myLim2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            testVal = matRes(i, m);
            Success = myComp(testVal, myLim) || myComp2(testVal, myLim2);
            
            if (Success)
                indexMatch.push_back(i);
        }
    }
    
    const int numCols = (keepRes) ? (m + 1) : m;
    const int numMatches = indexMatch.size();
    
    if (bLower)
        nRows = numMatches;
    else
        nRows  = (numMatches > userRows) ? userRows : numMatches;
    
    typeRcpp returnMatrix = Rcpp::no_init_matrix(nRows, numCols);
    const int lastCol = keepRes ? (m + 1) : m;
    
    for (int i = 0; i < nRows; ++i)
        for (int j = 0; j < lastCol; ++j)
            returnMatrix(i, j) = matRes(indexMatch[i], j);
    
    return returnMatrix;
}

template <typename typeVector>
void ApplyFunction(int n, int m, typeVector sexpVec, bool IsRep, int nRows, bool IsComb,
                   std::vector<int> myReps, SEXP ans, std::vector<int> freqs,
                   std::vector<int> z, bool IsMultiset, SEXP sexpFun, SEXP rho, int count) {
    if (IsComb) {
        if (IsMultiset)
            MultisetComboApplyFun(n, m, sexpVec, myReps, freqs, nRows, z, count, sexpFun, rho, ans);
        else
            ComboGeneralApplyFun(n , m, sexpVec, IsRep, count, nRows, z, sexpFun, rho, ans);
    } else {
        PermutationApplyFun(n, m, sexpVec, IsRep,nRows, IsMultiset, z, count, sexpFun, rho, ans);
    }
}

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool checkIsInteger(std::string funPass, unsigned long int uM, int n,
                    std::vector<double> rowVec, std::vector<double> vNum,
                    std::vector<double> myLim, funcPtr<double> myFunDbl,
                    bool checkLim = false) {
    
    if (funPass == "mean")
        return false;
    
    std::vector<double> vAbs;
    for (int i = 0; i < n; ++i)
        vAbs.push_back(std::abs(vNum[i]));
    
    double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    for (std::size_t i = 0; i < uM; ++i)
        rowVec[i] = static_cast<double>(vecMax);
    
    double testIfInt = myFunDbl(rowVec, uM);
    if (testIfInt > std::numeric_limits<int>::max())
        return false;
    
    if (checkLim) {
        vAbs.clear();
        for (std::size_t i = 0; i < myLim.size(); ++i)
            vAbs.push_back(std::abs(myLim[i]));
        
        double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
        if (vecMax > std::numeric_limits<int>::max())
            return false;
    }
    
    return true;
}

void getStartZ(int n, int m, double &lower, int stepSize, mpz_t &myIndex, bool IsRep,
               bool IsComb, bool IsMultiset, bool isGmp, std::vector<int> &myReps,
               std::vector<int> &freqsExpanded, std::vector<int> &startZ) {
    
    if (isGmp) {
        mpz_add_ui(myIndex, myIndex, stepSize);
        if (IsComb)
            startZ = nthCombinationGmp(n, m, myIndex, IsRep, IsMultiset, myReps);
        else
            startZ = nthPermutationGmp(n, m, myIndex, IsRep, IsMultiset, myReps, freqsExpanded, true);
    } else {
        lower += stepSize;
        if (IsComb)
            startZ = nthCombination(n, m, lower, IsRep, IsMultiset, myReps);
        else
            startZ = nthPermutation(n, m, lower, IsRep, IsMultiset, myReps, freqsExpanded, true);
    }
}

#endif
