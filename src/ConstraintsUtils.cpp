#include "Constraints/ConstraintsUtils.h"
#include "CleanConvert.h"

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool CheckIsInteger(const std::string &funPass, int n,
                    int m, const std::vector<double> &vNum,
                    const std::vector<double> &targetVals,
                    const funcPtr<double> myFunDbl, bool checkLim) {
    
    if (funPass == "mean") {
        return false;
    }
    
    std::vector<double> vAbs;
    
    for (int i = 0; i < n; ++i) {
        vAbs.push_back(std::abs(vNum[i]));
    }
    
    const double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    const std::vector<double> rowVec(m, vecMax);
    const double testIfInt = myFunDbl(rowVec, static_cast<std::size_t>(m));
    
    if (testIfInt > std::numeric_limits<int>::max()) {
        return false;
    }
    
    if (checkLim) {
        vAbs.clear();
        
        for (std::size_t i = 0; i < targetVals.size(); ++i) {
            if (static_cast<std::int64_t>(targetVals[i]) != targetVals[i]) {
                return false;
            } else {
                vAbs.push_back(std::abs(targetVals[i]));
            }
        }
        
        const double limMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
        
        if (limMax > std::numeric_limits<int>::max()) {
            return false;
        }
    }
    
    return true;
}

void fun(std::vector<double> &targetVals, std::vector<std::string> &compFunVec,
         SEXP Rtarget, SEXP RcompFun, bool IsConstrained) {
    
    // if (IsConstrained) {                // numOnly = true, checkWhole = false, negPoss = true
    //     CleanConvert::convertVector(Rtarget, targetVals, VecType::Numeric
    //                                 "limitConstraints",
    //                                 true, false, true);
    // 
    //     int len_comp = Rf_length(RcompFun);
    // 
    //     for ()
    // 
    // 
    //     s = malloc(len_comp + 1);
    //     if (s != NULL) {
    //         for (i = 0; i < n; i++) {
    //             s[i] = *CHAR(STRING_ELT(chars, i));
    //         }
    //         s[n] = '\0';
    //     } else {
    //         /*handle malloc failure*/
    //     }
    // 
    //     compFunVec = Rcpp::as<std::vector<std::string>>(f2);
    // 
    //     bool IsBetweenComp = false;
    //     ConstraintSetup(compFunVec, targetVals, IsBetweenComp);
    // 
    //     if (myType == VecType::Integer && !CheckIsInteger(mainFun, n, m, vNum, targetVals, funDbl, true)) {
    //         myType = VecType::Numeric;
    //     }
    // 
    //     double tolerance = 0;
    //     AdjustTargetVals(n, myType, targetVals, targetIntVals,
    //                      Rtolerance, compFunVec, tolerance, mainFun, vNum);
    // 
    //     if (myType == VecType::Integer) {
    //         GetPartitionCase(compFunVec, vInt, mainFun, targetIntVals, mySign,
    //                          PartType, ConstType, distinctTest, Rlow, myReps, n, m,
    //                          tolerance, IsMult, IsRep, IsBetweenComp, Rf_isNull(Rm));
    //     } else {
    //         GetPartitionCase(compFunVec, vNum, mainFun, targetVals, mySign,
    //                          PartType, ConstType, distinctTest, Rlow, myReps, n, m,
    //                          tolerance, IsMult, IsRep, IsBetweenComp, Rf_isNull(Rm));
    //     }
    // }
}
