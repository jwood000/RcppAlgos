#include "UserConstraintFuns.h"
#include "CombPermResultPtr.h"

// This is called when we can't easily produce a (loose) monotonic sequence overall,
// and we must generate and test every possible combination/permutation. This occurs
// when we are using "prod" and we have negative numbers involved. We also call this
// when lower is invoked implying that we are testing a specific range.
template <typename typeRcpp, typename T>
typeRcpp ConstraintsSpecial(int n, int m, std::vector<T> &v, bool IsRep, int nRows, 
                            bool keepRes, std::vector<int> z, double lower, const std::string &mainFun,
                            bool IsMult, double computedRows, const std::vector<std::string> &compFunVec,
                            const std::vector<T> &targetVals, bool IsComb, const std::vector<int> &freqs,
                            bool bLower, double userRows) {
    
    if (!bLower) {
        if (computedRows > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
        
        nRows = static_cast<int>(computedRows);
    }
    
    std::vector<int> indexMatch;
    indexMatch.reserve(nRows);
    
    Rcpp::XPtr<funcPtr<T>> xpFun = putFunPtrInXPtr<T>(mainFun);
    funcPtr<T> myFun = *xpFun;
    typeRcpp matRes = Rcpp::no_init_matrix(nRows, m + 1);
    Rcpp::XPtr<combPermResPtr<typeRcpp, T>> xpFunCoPeResPtr = 
        putCombResPtrInXPtr<typeRcpp, T>(IsComb, IsMult, IsRep);
    
    const combPermResPtr<typeRcpp, T> myFunCombPerm = *xpFunCoPeResPtr;
    myFunCombPerm(matRes, v, z, n, m, 0, nRows, freqs, myFun);
    
    Rcpp::XPtr<compPtr<T>> xpComp = putCompPtrInXPtr<T>(compFunVec[0]);
    compPtr<T> myComp = *xpComp;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            const T testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals))
                indexMatch.push_back(i);
        }
    } else {
        Rcpp::XPtr<compPtr<T>> xpComp2 = putCompPtrInXPtr<T>(compFunVec[1]);
        compPtr<T> myComp2 = *xpComp2;
        std::vector<T> targetVals2 = targetVals;
        targetVals2.erase(targetVals2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            const T testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals) || myComp2(testVal, targetVals2))
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

template Rcpp::IntegerMatrix ConstraintsSpecial(int, int, std::vector<int>&, bool, int, bool,
                                                std::vector<int>, double, const std::string&, bool, double,
                                                const std::vector<std::string>&, const std::vector<int>&, bool,
                                                const std::vector<int>&, bool, double);

template Rcpp::NumericMatrix ConstraintsSpecial(int, int, std::vector<double>&, bool, int, bool,
                                                std::vector<int>, double, const std::string&, bool, double,
                                                const std::vector<std::string>&, const std::vector<double>&, bool,
                                                const std::vector<int>&, bool, double);

