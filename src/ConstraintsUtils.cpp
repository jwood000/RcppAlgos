#include "ConstraintsUtils.h"
#include <cmath>

void BinaryNextElem(int &uppBnd, int &lowBnd, int &ind, int lastElem,
                    std::int64_t target, std::int64_t partial,
                    const std::vector<std::int64_t> &v) {
    
    std::int64_t dist = target - (partial + v[ind]);
    
    while ((uppBnd - lowBnd) > 1 && dist != 0) {
        const int mid = (uppBnd - lowBnd) / 2;
        ind = lowBnd + mid;
        dist = target - (partial + v[ind]);
        
        if (dist > 0)
            lowBnd = ind;
        else
            uppBnd = ind;
    }
    
    // Check last index. N.B. There are some cases when ind == lowBnd and dist < 0.
    // This will not matter as we simply reassign ind and recompute dist 
    if (dist < 0) {
        ind = lowBnd;
        dist = target - (partial + v[ind]);
    }
    
    // We must have dist <= 0. Below is an informal proof.
    // The sub-sequences are defined as below:
    //                  A_max = a_(i + 1), a_(i + 2), ..., a_m
    //                  A_set = a_1, a_2, ..., a_(i - 1)
    // A_set are those elements that have already been determined by the algorithm.
    // A_max is maximal (i.e. constructed of the the (m - i) largest elements). We
    // seek to determine the i_th element given the following contraints:
    //                      A_sum = A_set + a_i + A_max
    //                       dist = target - A_sum
    // With the goal of finding the minimum lexicographic combination such that the
    // dist = 0 (i.e. target = A_sum). If we have dist > 0 for any i, then it will
    // be impossible to obtain dist = 0. dist > 0 implies that the target > A_sum,
    // and since A_max is already maximal, we are not able to increase A_sum in
    // later iterations, thus we must have dist <= 0 for all i.
    if (dist > 0 && ind < lastElem)
        ++ind;
}

int GetFirstCombo(int m, const std::vector<std::int64_t> &v, bool IsRep, bool IsMult,
                  std::vector<int> &z, const std::vector<int> &freqs, std::int64_t target,
                  std::vector<int> repsCounter, int lastCol, int lastElem) {
    
    std::int64_t testMax = 0;
    constexpr std::int64_t zero64 = 0;
    
    if (IsRep) {
        testMax = v[lastElem] * m;
    } else if (IsMult) {
        const int lenMinusM = freqs.size() - m;
        
        for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j){
            testMax += v[freqs[i]];
        }
        
    } else {
        testMax = std::accumulate(v.cend() - m, v.cend(), zero64);
    }
    
    if (testMax < target)  {return 0;}
    int zExpCurrPos = IsMult ? freqs.size() - m : 0;
    int currPos = IsMult ? freqs[zExpCurrPos] : (IsRep ? lastElem : (v.size() - m));
    
    std::int64_t partial = testMax;
    partial -= v[currPos];
    std::int64_t testMin = 0;
    
    if (IsRep) {
        testMin = v[0] * m;
    } else if (IsMult) {
        
        for (int i = 0; i < m; ++i) {
            testMin += v[freqs[i]];
        }
        
    } else {
        testMin = std::accumulate(v.cbegin(), v.cbegin() + m, zero64);
    }
    
    if (testMin > target) {return 0;}
    int mid = currPos / 2;
    std::int64_t dist = target - (partial + v[mid]);
    
    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? currPos : mid;
    int ind = mid;
    
    for (int i = 0; i < lastCol; ++i) {
        BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
        z[i] = ind;
        partial += v[ind];
        
        if (IsMult) {
            --repsCounter[ind];
            
            if (repsCounter[ind] == 0)
                ++ind;
            
            ++zExpCurrPos;
            currPos = freqs[zExpCurrPos];
        } else if (!IsRep) {
            ++ind;
            ++currPos;
        }
        
        lowBnd = ind;
        uppBnd = currPos;
        mid = (uppBnd - lowBnd) / 2;
        
        ind = lowBnd + mid;
        partial -= v[currPos];
    }
    
    BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
    z[lastCol] = ind;
    
    // The algorithm above finds the first possible sum that equals
    // target. If there is no combination of elements from v that sum
    // to target, the algo returns the combination such that its sum
    // is closest to target and greater than target
    std::int64_t finalCheck = 0;
    
    for (int i = 0; i < m; ++i)
        finalCheck += v[z[i]];
    
    if (finalCheck != target)
        return 0;
    
    return 1;
}

template <typename typeVector>
inline void PopulateVec(int m, const std::vector<typeVector> &v,
                        std::vector<int> &z, int &count, int maxRows,
                        bool IsComb, std::vector<typeVector> &combinatoricsVec) {
    
    if (IsComb) {
        for (int k = 0; k < m; ++k)
            combinatoricsVec.push_back(v[z[k]]);
        
        ++count;
    } else {
        do {
            for (int k = 0; k < m; ++k)
                combinatoricsVec.push_back(v[z[k]]);
            
            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
    }
}

template <typename typeVector>
void SectionOne(const std::vector<typeVector> &v, std::vector<typeVector> &testVec,
                std::vector<int> &z, const std::vector<typeVector> &targetVals,
                std::vector<typeVector> &combinatoricsVec, std::vector<typeVector> &resultsVec,
                bool &t_0, bool &t_1, int &count, partialPtr<typeVector> partialFun,
                funcPtr<typeVector> constraintFun, compPtr<typeVector> compFunOne,
                compPtr<typeVector> compFunTwo, int m, int m1, int maxRows,
                int maxZ, bool IsComb, bool xtraCol) {
    
    for (int i = 0; i < m; ++i)
        testVec[i] = v[z[i]];
    
    const typeVector partialVal = constraintFun(testVec, m1);
    typeVector testVal = partialFun(partialVal, testVec.back(), m);
    t_0 = compFunTwo(testVal, targetVals);
    
    while (t_0 && t_1) {
        if (compFunOne(testVal, targetVals)) {
            int myStart = count;
            PopulateVec(m, v, z, count, maxRows, IsComb, combinatoricsVec);
            
            if (xtraCol)
                for (int i = myStart; i < count; ++i)
                    resultsVec.push_back(testVal);
            
            t_1 = count < maxRows;
        }
        
        t_0 = z[m1] != maxZ;
        
        if (t_0) {
            ++z[m1];
            testVec[m1] = v[z[m1]];
            testVal = partialFun(partialVal, testVec.back(), m);
            t_0 = compFunTwo(testVal, targetVals);
        }
    }
}

// This one is similar to DistinctAttr, however since v is non-standard,
// more work is required to obtain the optimal width. We cannot use the
// constant time formula and thus must resort to an O(n) algo.
distinctType DistinctAttrMapped(const std::vector<std::int64_t> &v, int m,
                                bool IsRep, bool IsMult, std::int64_t target, 
                                const std::vector<int> &Reps, Sign mySign,
                                bool IncludeZero, bool mIsNull) {
    int limit = 0;
    bool getAll = false;
    
    if (IsMult || !IsRep) {
        // No fancy tricks here... just sum until we breach
        std::int64_t testTar = 0;
        
        if (mySign == Sign::Positive) {
            for (std::size_t i = IncludeZero; i < v.size() && testTar <= target; ++i) {
                testTar += v[i];
                ++limit;
            }
        } else if (mySign == Sign::Negitive) {
            const int strt = (IncludeZero) ? v.size() - 2 : v.size() - 1;
            
            for (int i = strt; i >= 0 && testTar >= target; --i) {
                testTar += v[i];
                ++limit;
            }
        } else {
            if (target < 0) {
                for (int i = v.size() - 1; i >= 0 && testTar >= target; --i) {
                    // We don't count the length added by a value of zero
                    if (v[i]) {
                        testTar += v[i];
                        ++limit;
                    }
                }
            } else {
                for (std::size_t i = 0; i < v.size() && testTar <= target; ++i) {
                    // We don't count the length added by a value of zero
                    if (v[i]) {
                        testTar += v[i];
                        ++limit;
                    }
                }
            }
        }
        
        if (testTar > target)
            --limit;
        
        if (IsMult) {
            // Ensure all elements except the element representing zero
            // (if zero exists) have a frequency of 1.
            
            bool allOne = true;
            
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (v[i]) {
                    if (Reps[i] != 1) {
                        allOne = false;
                    }
                }
            }
            
            // The only way that getAll can be true in the non-traditional
            // way, is for our case to be a positive multiple of a the
            // traditional "getAll" scenario. E.g.
            // Tradtional "getAll":
            //     comboGeneral(0:10, freqs = c(10, rep(1, 10)),
            //                  constraintFun = "sum",
            //                  comparisonFun = "==",
            //                  limitConstraints = 10)
            //
            // Above example multiplied by 10:
            //     comboGeneral(seq(0L, 100L, 10L), m = 4,
            //                  freqs = c(10, rep(1, 10)),
            //                  constraintFun = "sum",
            //                  comparisonFun = "==",
            //                  limitConstraints = 100)
            //
            // If we included negative multiples, the result would not be
            // in lexicographical order.
            //
            // The result below is correct and clearly does not map to
            // the traditional case:
            //     comboGeneral(seq(0L, -100L, -10L), 4,
            //                  freqs = c(10, rep(1, 10)),
            //                  constraintFun = "sum",
            //                  comparisonFun = "==",
            //                  limitConstraints = -100)
            //          [,1] [,2] [,3] [,4]
            //     [1,] -100    0    0    0
            //     [2,]  -90  -10    0    0
            //     [3,]  -80  -20    0    0
            //     [4,]  -70  -30    0    0
            //     [5,]  -70  -20  -10    0
            //     [6,]  -60  -40    0    0
            //     [7,]  -60  -30  -10    0
            //     [8,]  -50  -40  -10    0
            //     [9,]  -50  -30  -20    0
            //    [10,]  -40  -30  -20  -10
            //
            // If we have all non-negative elements, IncludeZero = true,
            // and allOne = true, we need to find a fourth condition
            // analogous to (lenV == target + 1) in the traditional test.
            // We discovered above that this only gets mapped properly
            // if we have the traditional case times a positive number.
            // Thus, the first positive element will have to be that
            // multiple. That is, 1 will be mapped to M and M will be
            // our multiple. To obtain the mapped target, we simply
            // divide our target by M.
            
            if (
                   IncludeZero
                && mySign == Sign::Positive
                && lenV == (target / v[1]) + 1
                && allOne
              )
            {
                if (m >= limit) {
                    getAll = (Reps.front() >= (limit - 1)) ? true : false;
                } else {
                    limit = m;
                }
            } else {
                // N.B. In the calling function we have ensured that if the
                // freqs arg is invoked with all ones, we set IsMult to false.
                // This means that this should not happen
                limit = 0;
            }
        } else if (!IsRep) { // I.e. all elements have frequency = 1
            if (m < limit) {
                limit = m;
            } else if (!mIsNull) {
                limit = 0;
            }
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

distinctType DistinctAttr(int lenV, int m, bool IsRep, bool IsMult, std::int64_t target,
                          const std::vector<int> &Reps, bool IncludeZero, bool mIsNull) {
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
        // (a = 1, b = 1, c = -2)
        //
        //      x = (-1 + sqrt(1 + 4 * 1 * 2n)) / 2 * 1
        //
        // After solving for x, if sum(1:x) > target, we know
        // that the solution with the fewest number of elements
        // will contain x - 1 elements, hence std::floor (or 
        // just integer division).
        
        limit = (-1 + std::sqrt(1.0 + 8.0 * target)) / 2;
        
        if (IsMult) {
            // Ensure all elements except the first element are 1. The first
            // element should be zero and thus have a higher frequency in
            // order to test for partitions of different length.
            
            const bool allOne = std::all_of(Reps.cbegin() + 1, Reps.cend(), 
                                            [](int v_i) {return v_i == 1;});
            
            // We need lenV == target + 1 because we could have a case where
            // we have zero and allOne, but we are missing at least one
            // element in order to guarantee we generate all possible
            // partitions: E.g. v = 0, 1, 2, 3, 4; Reps = c(4, rep(1, 4));
            // and target = 6
            if (   
                   IncludeZero
                && lenV == target + 1
                && allOne
               )
            {
                if (m >= limit) {
                    getAll = (Reps.front() >= (limit - 1));
                } else {
                    limit = m;
                }
            } else {
                // N.B. In the calling function we have ensured that if the 
                // freqs arg is invoked with all ones, we set IsMult to false.
                // This means that this should not happen
                limit = 0;
            }
        } else if (!IsRep) {
            if (m < limit) {
                limit = m;
            } else if (!mIsNull) {
                limit = 0;
            }
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
            if (static_cast<std::int64_t>(targetVals[i]) != targetVals[i])
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
                        const std::vector<std::int64_t> &v, std::vector<int> &z,
                        const std::vector<int> &Reps, int target, int lenV,
                        int m, bool IncludeZero, bool IsRep, bool IsMult) {
    
    const bool IsDistinct = (PartType >= PartitionType::DstctStdAll);
    const int lastCol = (IsDistinct) ? distinctTest.limit - 1 : m - 1;
    const std::size_t partWidth = lastCol + 1;
    
    // Resize z appropriately
    z.assign(partWidth, 0);
    
    switch (PartType) {
        case PartitionType::Traditional: {
            z[lastCol] = target;
            break;
        } case PartitionType::TradNoZero: {
            std::fill(z.begin(), z.end(), 1);
            z[lastCol] = target - lastCol;
            break;
        } case PartitionType::DstctStdAll: {
            z[lastCol] = target;
            break;
        } case PartitionType::DstctShort: {
            z[lastCol] = target;
            break;
        } case PartitionType::DstctSpecial: {
            std::iota(z.begin() + Reps[0] - 1, z.end(), 0);
            z[lastCol] = target - (partWidth - Reps[0]) * (partWidth - Reps[0] - 1) / 2;
            break;
        } case PartitionType::DstctOneZero: {
            std::iota(z.begin(), z.end(), 0);
            z[lastCol] = target - (partWidth - 1) * (partWidth - 2) / 2;
            break;
        } case PartitionType::DstctNoZero: {
            std::iota(z.begin(), z.end(), 1);
            z[lastCol] = target - partWidth * (partWidth - 1) / 2;
            break;
        } default: {
            // Do something TBD
            break;
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
    
    if (compFunVec.size() > 2)
        Rcpp::stop("there cannot be more than 2 comparison operators");
    
    // The first 5 are "standard" whereas the 6th and 7th
    // are written with the equality first. Converting
    // them here makes it easier to deal with later.
    for (std::size_t i = 0; i < compFunVec.size(); ++i) {
        auto itComp = compForms.find(compFunVec[i]);
        
        if (itComp == compForms.end()) {
            Rcpp::stop("comparison operators must be one of the following: "
                           "'>', '>=', '<', '<=', or '=='");
        } else {
            compFunVec[i] = itComp->second;
        }
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
            // compSpecial in ConstraintsMain.h. If we look at the definitions of these comparison
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
                if (static_cast<std::int64_t>(vNum[i]) != vNum[i])
                    IsWhole = false;
                
            for (std::size_t i = 0; i < targetVals.size() && IsWhole; ++i)
                if (static_cast<std::int64_t>(targetVals[i]) != targetVals[i])
                    IsWhole = false;
                    
            tolerance = (IsWhole && mainFun != "mean") ? 0 : defaultTolerance;
        } else {
            // numOnly = true, checkWhole = false, negPoss = false, decimalFraction = true
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

template <typename typeVector>
int GetMappedTarget(const std::vector<typeVector> &v,
                    std::int64_t target, const std::vector<int> &Reps,
                    int m, int lenV, bool IsMult, bool IsRep) {
    
    std::vector<int> z(m, 0);
    std::vector<int> zExpanded;
    int myTarget = 0;
    
    std::vector<std::int64_t> v64(v.cbegin(), v.cend());
    const int lastCol = m - 1;
    
    for (std::size_t i = 0; i < Reps.size(); ++i)
        for (int j = 0; j < Reps[i]; ++j)
            zExpanded.push_back(i);
    
    Rcpp::Rcout << "\n";
    Rcpp::print(Rcpp::wrap(target));
    Rcpp::print(Rcpp::wrap(lastCol));
    Rcpp::print(Rcpp::wrap(lenV - 1));
    Rcpp::print(Rcpp::wrap(m));
    Rcpp::print(Rcpp::wrap(v64));
    
    if (GetFirstCombo(m, v64, IsRep, IsMult, z, zExpanded,
                      target, Reps, lastCol, lenV - 1)) {
        
        Rcpp::Rcout << "\n";
        Rcpp::print(Rcpp::wrap("First combo index"));
        Rcpp::print(Rcpp::wrap(z));
        
        myTarget = std::accumulate(z.begin(), z.end(), 0) + m;
    }
    
    return myTarget;
}

// IsBet is short for IsBetween
template <typename typeVector>
void GetPartitionCase(const std::vector<std::string> &compFunVec, std::vector<typeVector> &v,
                      const std::string &mainFun, const std::vector<typeVector> &target,
                      Sign mySign, PartitionType &PartType, ConstraintType &ConstType, Sign mySign,
                      distinctType &distinctTest, const SEXP &Rlow, std::vector<int> &Reps,
                      int lenV, int &m, double tolerance, bool IsMult,
                      bool IsRep, bool IsBet, bool mIsNull) {
    
    ConstType = ConstraintType::General;
    bool bLower = false;
    
    // Currently, we are not able to generate the nth
    // lexicographical partition. Thus, if lower is
    // non-trivial, we must use most general algo.
    if (!Rf_isNull(Rlow)) {
        auto tempLower = FromCpp14::make_unique<mpz_t[]>(1);
        mpz_init(tempLower[0]);
        
        createMPZArray(Rlow, tempLower.get(), 1, "lower");
        bLower = mpz_cmp_si(tempLower[0], 1) > 0;
    }
    
    // compFunVec should be non-empty if we made it this far.
    // Doesn't hurt to check
    if (!compFunVec.empty() && !bLower) {
        
        /// We start by assuming we don't have a nice partition case
        bool PartitionCase = false;
        
        if (IsMult) {
            for (int i = 0; i < (lenV - 1); ++i) {
                for (int j = (i + 1); j < lenV; ++j) {
                    if (v[i] > v[j]) {
                        std::swap(v[i], v[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(v.begin(), v.end());
        }
        
        std::int64_t tarTest = 0;
        
        if (   compFunVec[0] == "=="
            &&       mainFun == "sum"
            &&          lenV >  1
            &&             m >  1
           )
        {
            // We need to make sure not to include zero in the check below.
            // This is so because the zero can be used to obtain partitions
            // of differing lengths. Under normal circumcstances, this is
            // no problem because we simply have 0, 1, 2,..., however with
            // the capped cases (i.e. they don't start at 0 or 1, e.g. 3:14)
            // this case would be excluded because (3 - 0) != (4 - 3). Note,
            // We can only do this if all values have the same sign. When we
            // have mixed signs numbers, indexing breaks down. We are making
            // the case that we can't get all partititions of every length
            // when this occurs. First off, we have mapping issues. E.g. Let
            // v: -15, -9, -3, 0, 3, 9,..., 99 and a target of 93. For m = 3,
            // the first partition is 0, 0, -15, 9, 99 which maps to 
            // 4, 4, 1, 6, 21 giving a new target of 37. Now let m = 5
            // The first partition is -15, -15, -15, 39, 99 which maps to 
            // 1, 1, 1, 11, 21 for a mapped target of 35 (which is not 37!).
            //
            // Secondly, even if we could map properly for differing lengths
            // we would have issues with ordering (lexicographically).
            
            std::vector<typeVector> pTest;
            
            if (mySign != Sign::MixedBag) {
                for (auto val: v)
                    if (val != 0)
                        pTest.push_back(val);
            } else {
                pTest = v;
            }
            
            std::sort(pTest.begin(), pTest.end());
            const typeVector tarDiff = pTest[1] - pTest[0];
            
            if (static_cast<std::int64_t>(pTest[0]) == pTest[0]) {
                PartitionCase = true;
                
                for (std::size_t i = 1; i < pTest.size(); ++i) {
                    const typeVector testDiff = pTest[i] - pTest[i - 1];
                    
                    if (       std::abs(testDiff - tarDiff)  > tolerance ||
                        static_cast<std::int64_t>(pTest[i]) != pTest[i]) {
                        
                        PartitionCase = false;
                        break;
                    }
                }
            }
            
            if (target.size() == 1 || target.front() == target.back()) {
                tarTest = static_cast<std::int64_t>(target.front());
                
                if (PartitionCase)
                    PartitionCase = (tarTest == target.front());
            } else {
                PartitionCase = false;
            }
        }
        
        if (PartitionCase) {
            // Now that we know we have partitions, we need to determine
            // if the final result needs to be mapped. There are a couple
            // of ways this can happen.
            // 
            // 1. If the first element isn't zero or one.
            // 2. If the distance between elements is greater than 1.
            //
            // The vector vBase will take on the underlying base partition.
            //
            // Note, we have already ensured above that if we have
            // negative values and the differenece between every value
            // isn't the same, then we don't meet the partition scenario.
            
            std::vector<int> vBase(v.size());
            int mappedTarget = 0;
            const typeVector testDiff = v[1] - v[0];
            
            const bool condition_1 = (v.front() != 1 && v.front() != 0) || testDiff != 1;
            
            if (condition_1 && !Sign::MixedBag) {
                // Set our constraint type to indicate mapping
                // will be needed. See PartitionMain.cpp
                ConstType = ConstraintType::PartMapping;
                std::iota(vBase.begin(), vBase.end(), 0);
                mappedTarget = GetMappedTarget(v, tarTest, Reps,
                                               m, lenV, IsMult, IsRep);
            } else {
                ConstType = ConstraintType::PartStandard;
                mappedTarget = tarTest;
                
                for (int i = 0; i < lenV; ++i)
                    vBase[i] = v[i];
            }
            
            Rcpp::Rcout << "\n";
            Rcpp::print(Rcpp::wrap("mappedVector"));
            Rcpp::print(Rcpp::wrap(vBase));
            Rcpp::Rcout << "\n";
            Rcpp::print(Rcpp::wrap("mappedTarget"));
            Rcpp::print(Rcpp::wrap(mappedTarget));
            Rcpp::stop("digg");
            
            // We sorted v above to ensure that the last element is the maximum
            const int myMax = vBase.back();
            const int IncludeZero = (vBase.front() == 0);
            
            // Remember, lenV is the length of the vector v, so we could have a
            // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
            // cause a problem if we were to allow this. We have already ensured 
            // that the distance between each element is the same. This means for
            // the example we gave, we would have length(unique(diff(v))) > 1,
            // which means PartitionCase would be false, and thus the general
            // algorithm would be executed.
            //
            // We do have to ensure that the smallest element is non-negative, othe-
            // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & limit = 10
            // would pass as v = 0:9, m = 7, rep = TRUE, & limit = 9, --or--
            // v = 1:10, m = 7, rep = TRUE, & limit = 10 (Hence v.front() >= 0)
            
            if (                
                   myMax == mappedTarget
                && lenV + IncludeZero == mappedTarget
                && vBase.front() >= 0
               )
            {
                distinctTest = DistinctAttr(lenV, m, IsRep, IsMult, mappedTarget,
                                            Reps, static_cast<bool>(IncludeZero), mIsNull);
                
                if (distinctTest.limit > 0) {
                    m = distinctTest.limit;
                    
                    if (IncludeZero) {
                        if (IsMult) {
                            if (distinctTest.getAll) {
                                PartType = PartitionType::DstctStdAll;
                            } else if (Reps[0] >= (m - 1)) {
                                PartType = PartitionType::DstctShort;
                            } else {
                                PartType = PartitionType::DstctSpecial;
                            }
                        } else {
                            PartType = PartitionType::DstctOneZero;
                        }
                    } else {
                        PartType = PartitionType::DstctNoZero;
                    }
                } else if (IsRep) {
                    if (IncludeZero) {
                        PartType = PartitionType::Traditional;
                    } else {
                        PartType = PartitionType::TradNoZero;
                    }
                    
                    if (m >= lenV) { 
                        if (IncludeZero) {
                            m = lenV - 1;
                        } else {
                            ConstType = ConstraintType::General;
                        }
                    }
                }
            }  else {
                if (IsRep) {
                    PartType = PartitionType::TradCapped;
                } else if (!IsMult) {
                    PartType = PartitionType::DistCapped;
                }
            }
        } else if (
                        (compFunVec[0] == "==" || IsBet)
                     && lenV > 1
                     && m > 1
                     && mainFun != "max"
                     && mainFun != "min"
                  )
        {
                       
            // N.B. When we arrive here, the user must provide the width.
            // That is, m cannot be NULL
            ConstType = ConstraintType::PartitionEsque;
        }
    }
}

template <typename typeVector>
bool CheckSpecialCase(int n, bool bLower, const std::string &mainFun,
                      const std::vector<typeVector> &vNum) {
    
    bool result = false;
    
    // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
    // will simply return (upper - lower) # of combinations/permutations that meet the criteria
    if (bLower) {
        result = true;
    } else if (mainFun == "prod") {
        for (int i = 0; i < n; ++i) {
            if (vNum[i] < 0) {
                result = true;
                break;
            }
        }
    }
    
    return result;
}

template void SectionOne(const std::vector<int>&, std::vector<int> &testVec,
                         std::vector<int>&, const std::vector<int>&,
                         std::vector<int>&, std::vector<int>&,
                         bool&, bool&, int&, partialPtr<int>, funcPtr<int>,
                         compPtr<int>, compPtr<int>, int, int, int, int, bool, bool);
template void SectionOne(const std::vector<double>&, std::vector<double> &testVec,
                         std::vector<int>&, const std::vector<double>&,
                         std::vector<double>&, std::vector<double>&,
                         bool&, bool&, int&, partialPtr<double>, funcPtr<double>,
                         compPtr<double>, compPtr<double>, int, int, int, int, bool, bool);

template void GetPartitionCase(const std::vector<std::string>&, std::vector<int>&,
                               const std::string&, const std::vector<int>&, PartitionType&,
                               ConstraintType&, distinctType&, const SEXP&, std::vector<int>&,
                               int, int&, double, bool, bool, bool, bool);
template void GetPartitionCase(const std::vector<std::string>&, std::vector<double>&,
                               const std::string&, const std::vector<double>&, PartitionType&,
                               ConstraintType&, distinctType&, const SEXP&, std::vector<int>&,
                               int, int&, double, bool, bool, bool, bool);

template bool CheckSpecialCase(int, bool, const std::string&, const std::vector<int>&);
template bool CheckSpecialCase(int, bool, const std::string&, const std::vector<double>&);
