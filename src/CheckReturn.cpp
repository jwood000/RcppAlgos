#include "SetUpUtils.h"
#include "CheckReturn.h"

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
            Rf_error("constraintFun must be passed as a character");

        if (!Rf_isString(f2))
            Rf_error("comparisonFun must be passed as a character");
    }

    return result;
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

SEXP CheckReturn(SEXP Rv, SEXP f1,
                 SEXP f2, SEXP Rtarget,
                 SEXP RKeepRes, SEXP stdFun) {

    int res = returnType::constraintFun;

    if (Rf_isNull(f1)) {
        res = returnType::standard;
    } else {
        if (Rf_isFactor(Rv)) {
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
        const bool applyFun = !Rf_isNull(stdFun) && !Rf_isFactor(Rv);

        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rf_error("FUN must be a function!");

            res = returnType::anonymousFun;
        }
    }

    SEXP sexpRes = Rf_ScalarInteger(res);
    return sexpRes;
}
