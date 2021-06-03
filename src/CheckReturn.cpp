#include "Constraints/ConstraintsUtils.h"
#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"
#include "CheckReturn.h"
#include "SetUpUtils.h"
#include <numeric>

bool CheckConstrnd(SEXP RCnstrntFun, SEXP RCompFun, SEXP Rtarget) {
    // No need to check myType, as we have already done
    // so in CheckStdRet. Same goes for IsFactor.
    bool result = !Rf_isNull(RCnstrntFun) &&
                  !Rf_isNull(RCompFun) &&
                  !Rf_isNull(Rtarget);

    if (result) {
        if (!Rf_isString(RCnstrntFun))
            Rf_error("constraintFun must be passed as a character");

        if (!Rf_isString(RCompFun))
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

enum CheckReturnType : int {
    constraintFun = 0,
    standard = 1,
    anonymousFun = 2,
};

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
                        bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");

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
            if (!Rf_isFunction(stdFun))
                Rf_error("FUN must be a function!");

            res = CheckReturnType::anonymousFun;
        }
    }

    SEXP sexpRes = Rf_ScalarInteger(res);
    return sexpRes;
}

enum CheckPartitionType : int {
    NotPartition = 0,
    Partition = 1,
    PartitionEsque = 2,
};

SEXP CheckPartition(SEXP Rv, SEXP Rm, SEXP RmainFun, SEXP RcompFun,
                    SEXP Rlow, SEXP Rtarget, SEXP Rtolerance) {
    
    bool IsBetween = false;
    bool IsPartition = false;
    bool bLower = false;
    
    std::vector<double> v;
    std::vector<double> target;
    double tolerance;
    
    int res = CheckPartitionType::NotPartition;
    
    CleanConvert::convertVector(Rtarget, target,
                                VecType::Numeric,
                                "limitConstraints",
                                true, false, true);
    
    std::vector<std::string> compFunVec;
    int len_comp = Rf_length(RcompFun);
    
    for (int i = 0; i < len_comp; ++i) {
        const std::string temp(CHAR(STRING_ELT(RcompFun, i)));
        compFunVec.push_back(temp);
    }
    
    // Currently, we are not able to generate the nth
    // lexicographical partition. Thus, if lower is
    // non-trivial, we must use the most general algo.
    if (!Rf_isNull(Rlow)) {
        auto tempLower = FromCpp14::make_unique<mpz_t[]>(1);
        mpz_init(tempLower[0]);
        
        createMPZArray(Rlow, tempLower.get(), 1, "lower");
        bLower = mpz_cmp_si(tempLower[0], 1) > 0;
    }
    
    VecType myType;
    SetType(myType, Rv);
    
    if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
        Rf_error("contraintFun must be one of the following:"
                     " 'prod', 'sum', 'mean', 'max', or 'min'");
    }
    
    const std::string mainFun(CHAR(STRING_ELT(RmainFun, 0)));
    const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), mainFun);
    
    if (funIt == mainFunSet.end()) {
        Rf_error("contraintFun must be one of the following:"
                     " 'prod', 'sum', 'mean', 'max', or 'min'");
    }
    
    if (myType == VecType::Integer || myType == VecType::Numeric) {
        if (Rf_length(Rv) == 1) {
            int seqEnd = 0;
            
            // numOnly = true, checkWhole = true, negPoss = true
            CleanConvert::convertPrimitive(Rv, seqEnd, VecType::Integer,
                                           "If v is not a character"
                                           " and of length 1, it",
                                           true, true, true);
            
            std::pair<int, int> mnmx = std::minmax(1, seqEnd);
            const int n = mnmx.second - mnmx.first + 1;
            constexpr int maxVecSize = std::numeric_limits<int>::max() / 2;
            
            if (n < maxVecSize) {
                v.resize(n);
            } else {
                Rf_error("Not enough memory! The vector you have"
                             " requested is larger than %s",
                             std::to_string(maxVecSize).c_str());
            }
            
            std::iota(v.begin(), v.end(), mnmx.first);
        } else {
            v = CleanConvert::GetNumVec<double>(Rv);
        }
        
        GetTolerance(v, target, mainFun, Rtolerance, tolerance);
        ConstraintStructure(compFunVec, target, IsBetween);
        
        // compFunVec should be non-empty if we made it this far.
        // Doesn't hurt to check
        if (!compFunVec.empty() && !bLower && v.size() > 1) {
            if (compFunVec[0] == "==" && mainFun == "sum") {
                if (static_cast<std::int64_t>(v[0]) == v[0]) {
                    
                    IsPartition = true;
                    const double tarDiff = v[1] - v[0];
                    
                    for (std::size_t i = 1; i < v.size(); ++i) {
                        const double testDiff = v[i] - v[i - 1];
                        
                        if (std::abs(testDiff - tarDiff)  > tolerance ||
                            static_cast<std::int64_t>(v[i]) != v[i]) {
                            
                            IsPartition = false;
                            break;
                        }
                    }
                }
                
                if (IsPartition &&
                    (target.size() == 1 || target.front() == target.back()) &&
                    static_cast<std::int64_t>(target.front()) == target.front()) {
                    
                    res = CheckPartitionType::Partition;
                } else {
                    IsPartition = false;
                }
            }
            
            if (!IsPartition &&
                (compFunVec[0] == "==" || IsBetween) &&
                mainFun != "max" &&
                mainFun != "min" &&
                !Rf_isNull(Rm)) {
                
                // N.B. When we arrive here, the user must provide the width.
                // That is, m cannot be NULL
                res = CheckPartitionType::PartitionEsque;
            }
        }
    }
    
    SEXP sexpRes = Rf_ScalarInteger(res);
    return sexpRes;
}
