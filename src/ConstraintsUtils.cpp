#include "ConstraintsUtils.h"

distinctType DistinctAttr(int lenV, int m, bool IsRep, bool IsMult, int64_t target,
                          const std::vector<int> &Reps, bool IncludeZero) {
    int limit = 0;
    bool getAll = false;
    
    if (IsMult || !IsRep) {
        // The eqn below can be derived by taking note that the
        // smallest number of elements whose sum is at least 
        // the target will be comprised of the first x numbers.
        // That is, we need to solve for x such that:
        //
        //        sum(1:(x - 1)) <= target <= sum(1:x)
        //
        // These are triangle numbers which have the form:
        //
        //              sum(1:x) = x * (x + 1) / 2
        //
        // Given n = target, we have:
        //
        //    x * (x + 1) / 2 >= n  -->>  x^2 + x - 2n >= 0
        //
        // Finally, using the quadratic formula, we obtain:
        //
        //      x = (-1 + sqrt(1 + 4 * 1 * 2n)) / 2 * 1
        //
        // After solving for x, if sum(1:x) > target, we know
        // that the solution with the fewest number of elements
        // will contain x - 1 elements, hence std::floor.
        
        limit = static_cast<int>(std::floor((-1 + std::sqrt(1.0 + 8.0 * target)) / 2));
        
        if (IsMult) {
            // Ensure all elements except the first element are 1. The first
            // element should be zero and thus have a higher frequency in
            // order to test for partitions of different length.
            
            const bool allOne = std::all_of(Reps.cbegin() + 1, Reps.cend(), 
                                            [](int v_i) {return v_i == 1;});
            
            if (IncludeZero && lenV == target + 1 && allOne) {
                if (m >= limit) {
                    getAll = (Reps.front() >= (limit - 1)) ? true : false;
                } else {
                    limit = m;
                }
            } else {
                // N.B. In the calling function we have ensured that if the 
                // freqs arg is invoked with all ones, we set IsMult to false.
                limit = 0;
            }
        } else if (!IsRep) {
            limit = std::min(m, limit);
        }
    }
    
    // N.B. if limit = 0, this means we either have IsRep = true,
    // or we are not going to use the optimized algorithm. In this
    // case, we revert to the general algorithm.
    
    distinctType res;
    res.limit = limit;
    res.getAll = getAll;
    
    return res;
}

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool CheckIsInteger(const std::string &funPass, int n, int m,
                    const std::vector<double> &vNum, const std::vector<double> &targetVals,
                    funcPtr<double> myFunDbl, bool checkLim) {
    
    if (funPass == "mean")
        return false;
    
    std::vector<double> vAbs;
    
    for (int i = 0; i < n; ++i)
        vAbs.push_back(std::abs(vNum[i]));
    
    const double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    const std::vector<double> rowVec(m, vecMax);
    const double testIfInt = myFunDbl(rowVec, static_cast<std::size_t>(m));
    
    if (testIfInt > std::numeric_limits<int>::max())
        return false;
    
    if (checkLim) {
        vAbs.clear();
        
        for (std::size_t i = 0; i < targetVals.size(); ++i) {
            if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                return false;
            else
                vAbs.push_back(std::abs(targetVals[i]));
        }
        
        const double limMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
        
        if (limMax > std::numeric_limits<int>::max())
            return false;
    }
    
    return true;
}

void SetStartPartitionZ(PartitionType PartType, distinctType distinctTest,
                        std::vector<int> &z, const std::vector<int> &Reps,
                        int target, int lenV, int m, bool IncludeZero) {
    
    const bool IsDistinct = (PartType > PartitionType::PartTradNoZero);
    const int lastCol =  (IsDistinct) ? distinctTest.limit - 1 : m - 1;
    const std::size_t partWidth = lastCol + 1;
    z.assign(lastCol + 1 ,0);
    
    switch (PartType) {
    case PartTraditional: {
        z[lastCol] = target;
        break;
    }
    case PartTradNoZero: {
        std::fill(z.begin(), z.end(), 1);
        z[lastCol] = target - partWidth + 1;
        break;
    }
    case PartDstctStdAll: {
        z[lastCol] = target;
        break;
    }
    case PartDstctShort: {
        z[lastCol] = target;
        break;
    }
    case PartDstctSpecial: {
        std::iota(z.begin() + Reps[0] - 1, z.end(), 0);
        z[lastCol] = target - (partWidth - Reps[0]) * (partWidth - Reps[0] - 1) / 2;
        break;
    }
    case PartDstctOneZero: {
        std::iota(z.begin(), z.end(), 0);
        z[lastCol] = target - (partWidth - 1) * (partWidth - 2) / 2;
        break;
    }
    case PartDstctNoZero: {
        std::iota(z.begin(), z.end(), 1);
        z[lastCol] = target - partWidth * (partWidth - 1) / 2;
        break;
    }
    default: {
        z[lastCol] = target;
    }
    }
}

void ConstraintSetup(std::vector<std::string> &compFunVec,
                     std::vector<double> &targetVals, bool &IsBetweenComp) {
    
    if (targetVals.size() > 2) {
        Rcpp::stop("there cannot be more than 2 limitConstraints");
    } else if (targetVals.size() == 2 && targetVals[0] == targetVals[1]) {
        Rcpp::stop("The limitConstraints must be different");
    }
    
    std::vector<std::string>::const_iterator itComp;
    
    if (compFunVec.size() > 2)
        Rcpp::stop("there cannot be more than 2 comparison operators");
    
    for (std::size_t i = 0; i < compFunVec.size(); ++i) {
        itComp = std::find(compForms.cbegin(), compForms.cend(), compFunVec[i]);
        
        if (itComp == compForms.end()) {
            Rcpp::stop("comparison operators must be one of the following: "
                           "'>', '>=', '<', '<=', or '=='");
        }
        
        int myIndex = std::distance(compForms.cbegin(), itComp);
        
        // The first 5 are "standard" whereas the 6th and 7th
        // are written with the equality first. Converting
        // them here makes it easier to deal with later.
        if (myIndex > 4) {myIndex -= 3;}
        compFunVec[i] = compForms[myIndex];
    }
    
    if (compFunVec.size() == 2) {
        if (targetVals.size() == 1) {
            compFunVec.pop_back();
        } else {
            if (compFunVec[0] == "==" || compFunVec[1] == "==") {
                Rcpp::stop("If comparing against two limitConstraints, the "
                               "equality comparisonFun (i.e. '==') cannot be used. "
                               "Instead, use '>=' or '<='.");
            }
            
            if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1)) {
                Rcpp::stop("Cannot have two 'less than' comparisonFuns or two 'greater than' "
                               "comparisonFuns (E.g. c('<', '<=') is not allowed).");
            }
            
            // The two cases below are for when we are looking for all combs/perms such that when
            // myFun is applied, the result is between 2 values. These comparisons are defined in
            // compSpecial in ConstraintsMaster.h. If we look at the definitions of these comparison
            // functions in ConstraintsUtils.h, we see the following trend for all 4 cases:
            //
            // bool greaterLess(stdType x, const std::vector<stdType> &y) {return (x < y[0]) && (x > y[1]);}
            //
            // That is, we need the maixmum value in targetVals to be the first value and the second
            // value needs to be the minimum. At this point, the constraint algorithm will be
            // identical to when comparisonFun = "==" (i.e. we allow the algorithm to continue
            // while the results are less than (or equal to in cases where strict inequalities are
            // enforced) the target value and stop once the result exceeds that value).
            
            if (compFunVec[0].substr(0, 1) == ">" && 
                std::min(targetVals[0], targetVals[1]) == targetVals[0]) {
                
                compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                IsBetweenComp = true;
            } else if (compFunVec[0].substr(0, 1) == "<" && 
                std::max(targetVals[0], targetVals[1]) == targetVals[0]) {
                
                compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                IsBetweenComp = true;
            }
            
            if (IsBetweenComp) {
                compFunVec.pop_back();
                
                if (std::max(targetVals[0], targetVals[1]) == targetVals[1])
                    std::swap(targetVals[0], targetVals[1]);
            }
        }
    } else {
        if (targetVals.size() == 2)
            targetVals.pop_back();
    }
}

void AdjustTargetVals(int n, VecType myType, std::vector<double> &targetVals,
                      std::vector<int> &targetIntVals, const SEXP &Rtolerance,
                      std::vector<std::string> &compFunVec, double &tolerance,
                      const std::string &mainFun, const std::vector<double> &vNum) {
    
    if (myType == VecType::Integer) {
        targetIntVals.assign(targetVals.cbegin(), targetVals.cend());
    } else {
        // We first check if we are getting double precision. If so, for the
        // non-strict inequalities, we must alter the limit by some epsilon:
        //
        //                 x <= y   --->>>   x <= y + e
        //                 x >= y   --->>>   x >= y - e
        //
        // Equality is a bit tricky as we need to check whether the absolute
        // value of the difference is less than epsilon. As a result, we
        // can't alter the limit with one alteration. Observe:
        //
        //          x == y  --->>>  |x - y| <= e , which gives:
        //
        //                      -e <= x - y <= e
        //
        //                    1.     x >= y - e
        //                    2.     x <= y + e
        //
        // As a result, we must define a specialized equality check for double
        // precision. It is 'equalDbl' and can be found in ConstraintsUtils.h
        
        if (Rf_isNull(Rtolerance)) {
            bool IsWhole = true;
            
            for (int i = 0; i < n && IsWhole; ++i)
                if (static_cast<int64_t>(vNum[i]) != vNum[i])
                    IsWhole = false;
                
                for (std::size_t i = 0; i < targetVals.size() && IsWhole; ++i)
                    if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                        IsWhole = false;
                    
                    tolerance = (IsWhole && mainFun != "mean") ? 0 : defaultTolerance;
        } else {
            CleanConvert::convertPrimitive(Rtolerance, tolerance, 
                                           "tolerance", true, false, false, true);
        }
        
        const auto itComp = std::find(compSpecial.cbegin(),
                                      compSpecial.cend(), compFunVec[0]);
        
        if (compFunVec[0] == "==") {
            targetVals.push_back(targetVals[0] - tolerance);
            targetVals[0] += tolerance;
        } else if (itComp != compSpecial.end()) {
            targetVals[0] += tolerance;
            targetVals[1] -= tolerance;
        } else if (compFunVec[0] == "<=") {
            targetVals[0] += tolerance;
        } else if (compFunVec[0] == ">=") {
            targetVals[0] -= tolerance;
        }
        
        if (compFunVec.size() > 1) {
            if (compFunVec[1] == "<=") {
                targetVals[1] += tolerance;
            } else if (compFunVec[1] == ">=") {
                targetVals[1] -= tolerance;
            }
        }
    }
}
