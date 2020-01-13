#include "StandardUtils.h"

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

// [[Rcpp::export]]
int CheckReturn(SEXP Rv, SEXP f1, SEXP f2, SEXP Rtarget,
                bool IsFactor, SEXP RKeepRes, SEXP stdFun) {
    
    int res = 0;
    
    if (Rf_isNull(f1)) {
        res = 1;
    } else {
        if (IsFactor) {
            res = 1;
        } else {
            VecType myType = VecType::Integer;
            SetType(myType, Rv);
            
            if (myType > VecType::Logical) {
                res = 1;
            } else {
                if (!Rf_isNull(f2) && !Rf_isNull(Rtarget)) {
                    res = 0;  // This is a constrained output
                } else if (Rf_isNull(f2) && Rf_isNull(Rtarget)) {
                    // This is applying a constrained func only
                    if (Rf_isNull(RKeepRes)) {
                        res = 0;
                    } else {
                        bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
                        
                        if (keepRes) {
                            res = 0;
                        } else {
                            res = 1;
                        }
                    }
                } else {
                    // This means the input is non-sensible... std res = only
                    res = 1;
                }
            }
        }
    }
    
    if (res) {
        const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
        
        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            res = 2;
        }
    }
    
    return res;
}
