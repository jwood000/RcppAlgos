#include "SetUpUtils.h"

bool CheckConstrnd(SEXP RCnstrntFun, SEXP RCompFun, SEXP Rtarget) {
    // No need to check myType, as we have already done
    // so in CheckStdRet. Same goes for IsFactor.
    bool result = !Rf_isNull(RCnstrntFun) &&
                  !Rf_isNull(RCompFun) &&
                  !Rf_isNull(Rtarget);

    if (result) {
        if (!Rf_isString(RCnstrntFun)) {
            cpp11::stop("constraintFun must be passed as a character");
        }

        if (!Rf_isString(RCompFun)) {
            cpp11::stop("comparisonFun must be passed as a character");
        }
    }

    return result;
}

[[cpp11::register]]
SEXP CheckConstrndCpp(SEXP RCnstrntFun, SEXP RCompFun, SEXP Rtarget) {
    return Rf_ScalarLogical(CheckConstrnd(RCnstrntFun, RCompFun, Rtarget));
}

// CheckReturn(v, constraintFun, comparisonFun,
//             limitConstraints, keepResults, FUN)
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
//      keepResults: Are we keeping the results of constraintFun?
//              FUN: This is an anonymous function

enum CheckReturnType : int {
    constraintFun = 0,
    standard = 1,
    anonymousFun = 2,
};

[[cpp11::register]]
SEXP CheckReturn(SEXP Rv, SEXP RCnstrntFun, SEXP RCompFun,
                 SEXP Rtarget, SEXP RKeepRes, SEXP stdFun) {

    int res = CheckReturnType::constraintFun;

    if (Rf_isNull(RCnstrntFun)) {
        res = CheckReturnType::standard;
    } else {
        if (Rf_isFactor(Rv)) {
            res = CheckReturnType::standard;
        } else {
            VecType myType = VecType::Integer;
            SetType(myType, Rv);

            if (myType > VecType::Numeric) {
                res = CheckReturnType::standard;
            } else {
                if (!Rf_isNull(RCompFun) && !Rf_isNull(Rtarget)) {
                    // This is a constrained output
                    res = CheckReturnType::constraintFun;
                } else if (Rf_isNull(RCompFun) && Rf_isNull(Rtarget)) {
                    // This is applying a constrained func only
                    if (Rf_isNull(RKeepRes)) {
                        res = CheckReturnType::constraintFun;
                    } else {
                        bool keepRes = CleanConvert::convertFlag(RKeepRes, "keepResults");

                        if (keepRes) {
                            res = CheckReturnType::constraintFun;
                        } else {
                            res = CheckReturnType::standard;
                        }
                    }
                } else {
                    // This means the input is non-sensible... std res = only
                    res = CheckReturnType::standard;
                }
            }
        }
    }

    // if res isn't constrained (i.e. CheckReturnType::constraintFun)
    if (res) {
        const bool applyFun = !Rf_isNull(stdFun) && !Rf_isFactor(Rv);

        if (applyFun) {
            if (!Rf_isFunction(stdFun)) {
                cpp11::stop("FUN must be a function!");
            }

            res = CheckReturnType::anonymousFun;
        }
    }

    return Rf_ScalarInteger(res);
}
