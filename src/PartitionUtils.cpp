#include "Partitions/PartitionsUtils.h"

void BinaryNextElem(int &uppBnd, int &lowBnd, int &ind, int lastElem,
                    std::int64_t target, std::int64_t partial,
                    const std::vector<std::int64_t> &v) {

    std::int64_t dist = target - (partial + v[ind]);

    while ((uppBnd - lowBnd) > 1 && dist != 0) {
        const int mid = (uppBnd - lowBnd) / 2;
        ind = lowBnd + mid;
        dist = target - (partial + v[ind]);

        if (dist > 0) {
            lowBnd = ind;
        } else {
            uppBnd = ind;
        }
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

    if (dist > 0 && ind < lastElem) {
        ++ind;
    }
}

// A return of -2 means that the length was too short.
// A return of -1 means the lengths is too long.
// A return of 0 means that the length was okay, however
//    no solution exists
// A return of 1 means we found a solution
int GetFirstPartition(const std::vector<std::int64_t> &v,
                      const std::vector<int> &freqs, std::vector<int> &z,
                      std::vector<int> repsCounter, std::int64_t target,
                      int m, int lastCol, int lastElem,
                      bool IsRep, bool IsMult) {

    std::int64_t testMax = 0;
    constexpr std::int64_t zero64 = 0;

    if (IsRep) {
        testMax = v[lastElem] * m;
    } else if (IsMult) {
        const int lenMinusM = freqs.size() - m;

        for (int i = freqs.size() - 1; i >= lenMinusM; --i){
            testMax += v[freqs[i]];
        }

    } else {
        testMax = std::accumulate(v.cend() - m, v.cend(), zero64);
    }
    
    // The length is too small
    if (testMax < target) {
        return -2;
    }
    
    int zExpCurrPos = IsMult ? freqs.size() - m : 0;
    int currPos = IsMult ? freqs[zExpCurrPos] : (IsRep ? lastElem : (v.size() - m));

    std::int64_t partial = testMax;
    partial -= v[currPos];
    std::int64_t testMin = 0;

    if (IsRep) {
        testMin = v.front() * m;
    } else if (IsMult) {

        for (int i = 0; i < m; ++i) {
            testMin += v[freqs[i]];
        }

    } else {
        testMin = std::accumulate(v.cbegin(), v.cbegin() + m, zero64);
    }
    
    // The length is too long
    if (testMin > target) {
        return -1;
    }
    
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

// ****************************************************************
// Let v_i = (v_i1, v_i2, ..., v_im) be the ith partition of our
// target, T. Thus v_i1 + v_i2 + ... + v_im = T. We want to find
// the mapped target value. For our case, the only mapped partitions
// occur by linear shifts (a.k.a. translations) and/or constant
// multiples. Thus we can forumlate a mapping function as:
//
//                   f(x) = (1 / a) * x - b
//
// Where a is the distance between consecutive non-zero elements in
// the original sorted vector and b is the distance of the smallest
// element in v_mapped from zero. This can be distilled to:
// 
// If IncludeZero = TRUE:
//
//               b = (smallest element of v) / a
//
// Else:
//
//             b = (smallest element of v) / a - 1
//             b = ((smallest element of v) - a) / a
//
// If we are not including zero, then we need to subtract 1.
// How we determine to include zero or not remains a mystery
// for now. I'm inclinded to always include zero as the mapping
// will be much easier. That is, we can index without constantly
// subtracting one.
//
// Note, this is true for positive, negative, and mixed cases.
//
// If we map each element of v_i, we obtain:
//
// v_i_mapped = ((v_i1 / a - b), (v_i2 / a - b), ..., (v_im / a - b))
//
// Thus our new mapped target T_mapped is:
//
// T_mapped = (v_i1 / a - b) + (v_i2 / a - b) + ... + (v_im / a - b)
// T_mapped = (v_i1 / a) + (v_i2 / a) + ... + (v_im / a) - b * m
// T_mapped = (1 / a) * (v_i1 + v_i2 + ... + v_im) - b * m
//
// This gives us our main result:
//
//                  T_mapped = (1 / a) * T - b * m
//
// Since we cannot guarantee that a will evenly divide T or the smalle
// element of v, we instead use the following form:
//
//                  T_mapped = (T - min(v) * m) / a
//
// This is true regardless of IsRep, IsMult, etc. as our original
// assumption is that we are given the ith partition.
// ****************************************************************

// This is similar to GetDesign, however since v is non-standard,
// more work is required to obtain the optimal width. We cannot use the
// constant time formula and thus must resort to an O(n) algo. The
// same assumptions that were stated in GetDesign hold here.
void GetDesignMapped(const std::vector<std::int64_t> &v,
                     const std::vector<int> &Reps, 
                     PartDesignType &part, int m,
                     int lenV, bool mIsNull) {
    
    int width = m;
    
    // We are guaranteed that v has at least 2 non-zero values.
    // Also note that v and Reps have already been sorted.
    if (part.sign == Sign::Positive) {
        if (v.front() == 0) {
            part.slope = v[2] - v[1];
        } else {
            part.slope = v[1] - v[0];
        }
    } else {
        // If we have Sign::Negative or Sign::MixedBag, the first element
        // is guaranteed to be non-zero. For Sign::Negative, since we know
        // there are at least 2 non-zero values, the second element must
        // also be non-zero. For Sign::MixedBag, all elements must be 
        // spaced equally regardless of the presence of zero or not.
        part.slope = v[1] - v[0];
    }
    
    if (part.isRep) {
        if (mIsNull) {
            part.getAll = true;
            part.includeZero = true;
            
            // If we end up here, we assume that we have zero since
            // we can find a mapping from 1:n; tar = t to 0:(n -1);
            // tar = t'. In this case our "smallest" element will be
            // the second element of v.
            width = std::abs(part.tar) / v[1];
        }
        
        part.mappedTar = (part.tar - part.shift * width) / part.slope;
    } else {
        // We have multisets or distinct case
        
        // Check if all elements except the first element are 1. For the
        // mapped case, the first element will probably not be zero,
        // however it will be isomorphic to zero if getAll = true.
        const bool allOne = (!part.isMult) ? false : 
                            std::all_of(Reps.cbegin() + 1, Reps.cend(),
                                        [](int v_i) {return v_i == 1;});
        
        part.includeZero = (Reps.front() > 1 && allOne);
        
        // Now, we sum until until we breach. How we accumulate depends
        // on the sign of target and the value of mySign
        std::int64_t testTar = 0;
        
        if (part.sign == Sign::Positive) {
            for (std::size_t i = part.includeZero; i < v.size() && testTar < part.tar; ++i) {
                for (int j = 0; j < Reps[i] && testTar < part.tar; ++j, ++width) {
                    testTar += v[i];
                }
            }
        } else if (part.sign == Sign::Negative) {
            for (std::size_t i = part.includeZero; i < v.size() && testTar > part.tar; ++i) {
                for (int j = 0; j < Reps[i] && testTar > part.tar; ++j, ++width) {
                    testTar += v[i];
                }
            }
        } else {
            // Again, we need to be careful here. If we have mixed signs, zero
            // will play a non-trivial role as it will eventually be mapped
            // to a non-zero value, so we must count it! Initially, we were
            // checking if (v[i] != 0) then count it. This was wrong! 
            //
            // For more info, check this example:
            //      comboGeneral(seq(-6, 24, 3), m = 4,
            //                   freqs = c(9, rep(1, 10)),
            //                   constraintFun = "sum",
            //                   comparisonFun = "==",
            //                   limitConstraints = 6)
            
            if (part.tar < 0) {
                for (std::size_t i = part.includeZero; i < v.size() && testTar > part.tar; ++i) {
                    for (int j = 0; j < Reps[i] && testTar > part.tar; ++j, ++width) {
                        testTar += v[i];
                    }
                }
            } else {
                for (std::size_t i = part.includeZero; i < v.size() && testTar < part.tar; ++i) {
                    for (int j = 0; j < Reps[i] && testTar < part.tar; ++j, ++width) {
                        testTar += v[i];
                    }
                }
            }
        }
        
        if (part.tar >= 0 && testTar > part.tar)
            --width;
        
        if (part.tar < 0 && testTar < part.tar)
            --width;
        
        // We must be careful when we are determining getAll
        //
        // Tradtional "getAll":
        //
        //     comboGeneral(0:10, freqs = c(10, rep(1, 10)),
        //                  constraintFun = "sum",
        //                  comparisonFun = "==",
        //                  limitConstraints = 10)
        //
        // traditional example multiplied by 10:
        //
        //     comboGeneral(seq(0L, 100L, 10L), m = 4,
        //                  freqs = c(10, rep(1, 10)),
        //                  constraintFun = "sum",
        //                  comparisonFun = "==",
        //                  limitConstraints = 100)
        //
        // traditional example shifted by 2:
        //
        //     comboGeneral(2:12, m = 4, 
        //                  freqs = c(10, rep(1, 10)),
        //                  constraintFun = "sum",
        //                  comparisonFun = "==",
        //                  limitConstraints = 18)
        //
        // traditional example shifted by 2 and multiplied by 10:
        //
        //     comboGeneral(seq(20, 120, 10), m = 4, 
        //                  freqs = c(10, rep(1, 10)),
        //                  constraintFun = "sum",
        //                  comparisonFun = "==",
        //                  limitConstraints = 180)
        //
        // 
        // We cannot be confused by the negative counterparts. That
        // is, -10:0, target = -10. This will not map appropriately
        // to the traditional case as the result must be in order.
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
        
        part.mappedTar = (part.tar - part.shift * width) / part.slope;
        
        // We need lenV == mappedTarget + 1 because we could have a case where
        // we have zero and allOne, but we are missing at least one
        // element in order to guarantee we generate all possible
        // partitions: E.g. v = 0, 1, 2, 3, 4; Reps = c(4, rep(1, 4));
        // and target = 6. This is equivalent to checking v.back() == target
        if (part.includeZero && lenV == part.mappedTar + 1 && allOne) {
            // We need at least 1 non-zero value
            // in order to set getAll = true. E.g
            //
            // Case: Distinct partitions of 10
            // width = 4 (determined above e.g. 1 2 3 4)
            //
            // 1st part:   0 0 0 10
            //
            // Thus zero must have at least a
            // frequency of 3 = 4 - 1
            
            if (mIsNull || m == width) {
                part.getAll = (Reps.front() >= (width - 1));
            }
        } else {
            part.mappedTar = (part.tar - part.shift * m) / part.slope;
        }
    }
    
    part.width = width;
}

// Our philosophy is that if the user supplies the width, we
// give that to them. If the user leaves it NULL, we figure
// out the 'best' width possible. You will see that we
// give preference to the longest possible width.
void GetDesign(const std::vector<int> &Reps, PartDesignType &part,
               int m, int lenV, bool mIsNull) {

    // Only holds when elements are distinct (except zero)
    //
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
    int width = (mIsNull) ? (-1 + std::sqrt(1.0 + 8.0 * part.tar)) / 2 : m;

    if (part.isMult) {
        // Ensure all elements except the first element are 1. The first
        // element should be zero and thus have a higher frequency in
        // order to test for partitions of different length.

        const bool allOne = std::all_of(Reps.cbegin() + 1, Reps.cend(),
                                        [](int v_i) {return v_i == 1;});

        // We need lenV == target + 1 because we could have a case where
        // we have zero and allOne, but we are missing at least one
        // element in order to guarantee we generate all possible
        // partitions: E.g. v = 0, 1, 2, 3, 4; Reps = c(4, rep(1, 4));
        // target = 6. This is equivalent to checking v.back() == target
        if (part.includeZero && lenV == part.tar + 1 && allOne) {
            // We need at least 1 non-zero value
            // in order to set getAll = true. E.g
            //
            // Case: Distinct partitions of 10
            // width = 4 (determined above)
            //
            // 1st part:   0 0 0 10
            //
            // Thus zero must have at least a
            // frequency of 3 = 4 - 1

            if (mIsNull || m == width) {
                part.getAll = (Reps.front() >= (width - 1));
            }
        } else {
            // N.B. In the calling function we have ensured that if the
            // freqs arg is invoked with all ones, we set IsMult to false.
            // From this point, we have standard multiset partition case

            if (mIsNull) {
                width = 0;

                // If IncludeZero = true, we start at index 1
                std::int64_t testTar = 0;

                for (int i = part.includeZero, val = 1;
                     i < lenV && testTar < part.tar; ++i, ++val) {

                    for (int j = 0; j < Reps[i] && testTar < part.tar; ++j, ++width) {
                        testTar += val;
                    }
                }

                if (testTar > part.tar)
                    --width;
            }
        }
    } else if (part.isRep && mIsNull) {
        width = part.tar; // i.e. 1 * target = target
    } else {
        // Do nothing... we've already found width above
    }

    part.width = width;
}

template <typename T>
int GetMappedTarget(const std::vector<T> &v, std::vector<int> &vBase,
                    const std::vector<int> &Reps, PartDesignType &part,
                    int m, int lenV) {

    std::vector<int> z(m, 0);
    std::vector<int> zExpanded;
    int myTarget = 0;

    std::vector<std::int64_t> v64(v.cbegin(), v.cend());
    const int lastCol = m - 1;

    for (std::size_t i = 0; i < Reps.size(); ++i)
        for (int j = 0; j < Reps[i]; ++j)
            zExpanded.push_back(i);
    
    const int res = GetFirstPartition(v64, zExpanded, z, Reps,
                                      part.tar, m, lastCol, lenV - 1,
                                      part.isRep, part.isMult);
    
    if (res == 1) {
        myTarget = std::accumulate(z.begin(), z.end(), 0) + m;
    }

    return myTarget;
}

template <typename T>
bool CheckPartition(const std::vector<std::string> &compFunVec,
                    const std::vector<T> &v, const std::string &mainFun,
                    const std::vector<T> &target, PartDesignType &part,
                    SEXP Rlow, int lenV, int m, double tolerance,
                    bool IsBetween, bool mIsNull) {

    /// We start by assuming we don't have a nice partition case
    part.ctype = ConstraintType::General;
    bool IsPartition = false;
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
    
    // only proceed if we have more than 1 non-zero value
    int numNonTrivElem = 0;
    
    for (auto v_i: v)
        if (v_i)
            ++numNonTrivElem;
    
    // compFunVec should be non-empty if we made it this far.
    // Doesn't hurt to check
    if (!compFunVec.empty() && !bLower && numNonTrivElem > 1 && m > 1) {
        
        std::int64_t tarTest = 0;

        if (compFunVec[0] == "==" && mainFun == "sum") {
            // We need to make sure not to include zero in the check below
            // when part.sign != Sign::MixedBag.
            //
            // This is so because the zero can be used to obtain partitions
            // of differing lengths. Under normal circumcstances, this is
            // no problem because we simply have 0, 1, 2,..., however with
            // the capped cases (i.e. they don't start at 0 or 1, e.g. 3:14)
            // this case would be excluded because (3 - 0) != (4 - 3). When we
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

            std::vector<T> pTest;

            if (part.sign != Sign::MixedBag) {
                for (auto val: v) {
                    if (val != 0) {
                        pTest.push_back(val);
                    }
                }
            } else {
                pTest = v;
            }

            std::sort(pTest.begin(), pTest.end());
            const T tarDiff = pTest[1] - pTest[0];

            if (static_cast<std::int64_t>(pTest[0]) == pTest[0]) {
                IsPartition = true;

                for (std::size_t i = 1; i < pTest.size(); ++i) {
                    const T testDiff = pTest[i] - pTest[i - 1];

                    if (std::abs(testDiff - tarDiff)  > tolerance ||
                        static_cast<std::int64_t>(pTest[i]) != pTest[i]) {

                        IsPartition = false;
                        break;
                    }
                }
            }

            if (target.size() == 1 || target.front() == target.back()) {
                tarTest = static_cast<std::int64_t>(target.front());

                if (IsPartition){
                    IsPartition = (tarTest == target.front());
                }
                
            } else {
                IsPartition = false;
            }
            
            if (IsPartition) {
                part.tar = tarTest;
            }
        }

        if (!IsPartition &&
            (compFunVec[0] == "==" || IsBetween) &&
            mainFun != "max" &&
            mainFun != "min" &&
            !mIsNull) {

            // N.B. When we arrive here, the user must provide the width.
            // That is, m cannot be NULL
            part.ctype = ConstraintType::PartitionEsque;
        }
    }

    return IsPartition;
}

#include <iostream>

template <typename T>
void GetPartitionCase(const std::vector<int> &Reps,
                      std::vector<T> &v, PartDesignType &part,
                      int lenV, int &m, bool mIsNull) {

    // Now that we know we have partitions, we need to determine
    // if we are in a mapping case. There are a couple of ways 
    // this can happen.
    //
    // 1. If the first element isn't zero or one.
    // 2. If the distance between elements is greater than 1.
    //
    // For point 2, we need to check the difference between two
    // elements of v. If the first element is 0, we have satisfied
    // the first point, and if the next 10 elements are
    // seq(1, 99, 10), then only checking the difference between 
    // the first 2 elements would not reveal the nonstandard
    // distance between each element.
    //
    // The vector vBase will take on the underlying base partition.
    //
    // Note, we have already ensured above that if we have
    // negative values and the differenece between every value
    // isn't the same, then we don't meet the partition scenario.

    std::vector<int> vBase(lenV);
    const T testDiff1 = v[1] - v[0];
    const T testDiff2 = v.size() > 2 ? v[2] - v[1] : 1;
    
    part.mappedTar = part.tar;
    const bool condition_1 = (v.front() != 1 && v.front() != 0) ||
                              testDiff1 != 1 || testDiff2 != 1;
    
    if (condition_1) {
        // Set our constraint type to indicate mapping will be
        // needed to determine # of partitions. See PartitionMain.cpp
        part.ctype = ConstraintType::PartMapping;
        
        // For right now, we require the m be provided by the
        // user. in the future our goal will be to determine
        // the width.
        
        if (mIsNull) {
            Rf_error("When finding non-canonical partitions,"
                         " you must provide the width, m.");
        }
        
        GetMappedTarget(v, vBase, Reps, part, m, lenV);
    } else {
        part.ctype = ConstraintType::PartStandard;

        for (int i = 0; i < lenV; ++i)
            vBase[i] = v[i];
    }
    
    std::cout << "mappedVector\n";
    
    for (auto v_i: vBase)
        std::cout << v_i << " ";

    std::cout << std::endl;
    
    std::cout << "\n mappedTarget " << part.mappedTar << std::endl;
    Rf_error("digg");

    // We sorted v above to ensure that the last element is the maximum
    const int myMax = vBase.back();

    // Remember, lenV is the length of the vector v, so we could have a
    // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
    // cause a problem if we were to allow this. We have already ensured
    // that the distance between each element is the same. This means for
    // the example we gave, we would have length(unique(diff(v))) > 1,
    // which means PartitionCase would be false, and thus the general
    // algorithm would be executed.
    //
    // We do have to ensure that the smallest element is non-negative, othe-
    // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & width = 10
    // would pass as v = 0:9, m = 7, rep = TRUE, & width = 9, --or--
    // v = 1:10, m = 7, rep = TRUE, & width = 10 (Hence v.front() >= 0)

    if (myMax == part.mappedTar &&
        lenV + part.includeZero == part.mappedTar &&
        vBase.front() >= 0) {
        
        // distinctTest = DistinctAttr(lenV, m, IsRep, IsMult, mappedTarget,
        //                             Reps, static_cast<bool>(IncludeZero), mIsNull);

        // if (distinctTest.width > 0) {
        //     m = distinctTest.width;
        // 
        //     if (IncludeZero) {
        //         if (IsMult) {
        //             if (distinctTest.getAll) {
        //                 PartType = PartitionType::DstctStdAll;
        //             } else if (Reps[0] >= (m - 1)) {
        //                 PartType = PartitionType::DstctShort;
        //             } else {
        //                 PartType = PartitionType::DstctSpecial;
        //             }
        //         } else {
        //             PartType = PartitionType::DstctOneZero;
        //         }
        //     } else {
        //         PartType = PartitionType::DstctNoZero;
        //     }
        // } else if (IsRep) {
        //     if (IncludeZero) {
        //         PartType = PartitionType::Traditional;
        //     } else {
        //         PartType = PartitionType::TradNoZero;
        //     }
        // 
        //     if (m >= lenV) {
        //         if (IncludeZero) {
        //             m = lenV - 1;
        //         } else {
        //             ConstType = ConstraintType::General;
        //         }
        //     }
        // }
    }  else {
        // if (IsRep) {
        //     PartType = PartitionType::TradCapped;
        // } else if (!IsMult) {
        //     PartType = PartitionType::DistCapped;
        // }
    }
}
