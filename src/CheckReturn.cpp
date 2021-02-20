#include "StandardUtils.h"

enum returnType : int {
    constraintFun = 0,
         standard = 1,
     anonymousFun = 2,
};

bool CheckConstrnd(SEXP f1, SEXP f2, SEXP Rtarget) {
    // No need to check myType, as we have already done 
    // so in CheckStdRet. Same goes for IsFactor. 
    bool result = !Rf_isNull(f1) && !Rf_isNull(f2) && !Rf_isNull(Rtarget);
    
    if (result) {
        if (!Rf_isString(f1))
            Rcpp::stop("constraintFun must be passed as a character");
        
        if (!Rf_isString(f2))
            Rcpp::stop("comparisonFun must be passed as a character");
    }
    
    return result;
}

// CheckReturn(v, constraintFun, comparisonFun,
//             limitConstraints, IsFactor, keepResults, FUN)
//
//                v: the source vector
//    constraintFun: "sum", "mean", "prod", "min", "max"
//    comparisonFun: "==", ">", "<", ">=", "<=" (Single comparisons)
//
// When using two comparison operations, we could have the following
// assuming limitConstraints also has two values with the first value
// less than the second value:
//
// **** Between: c(">", "<" ), c(">", "<="), c(">=", "<"), c(">=", "<=")
// **** Outside: c("<",  ">"), c("<=", ">=")
//
// limitConstraints: The value(s) that will be used for comparison.
//         IsFactor: Are we dealing with factors?
//      keepResults: Are we keeping the results of constraintFun?
//              FUN: This is an anonymous function 

// [[Rcpp::export]]
int CheckReturn(SEXP Rv, SEXP f1, SEXP f2, SEXP Rtarget,
                bool IsFactor, SEXP RKeepRes, SEXP stdFun) {
    
    int res = returnType::constraintFun;
    
    if (Rf_isNull(f1)) {
        res = returnType::standard;
    } else {
        if (IsFactor) {
            res = returnType::standard;
        } else {
            VecType myType = VecType::Integer;
            SetType(myType, Rv);
            
            if (myType > VecType::Numeric) {
                res = returnType::standard;
            } else {
                if (!Rf_isNull(f2) && !Rf_isNull(Rtarget)) {
                    res = returnType::constraintFun;  // This is a constrained output
                } else if (Rf_isNull(f2) && Rf_isNull(Rtarget)) {
                    // This is applying a constrained func only
                    if (Rf_isNull(RKeepRes)) {
                        res = returnType::constraintFun;
                    } else {
                        bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
                        
                        if (keepRes) {
                            res = returnType::constraintFun;
                        } else {
                            res = returnType::standard;
                        }
                    }
                } else {
                    // This means the input is non-sensible... std res = only
                    res = returnType::standard;
                }
            }
        }
    }
    
    // if res isn't constrained (i.e. returnType::constraintFun)
    if (res) {
        const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
        
        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            res = returnType::anonymousFun;
        }
    }
    
    return res;
}
