#include "cpp11/protect.hpp"
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
                      int m, int lastCol, int lenV, bool IsRep, bool IsMult) {

    std::int64_t testMax = 0;
    const int lastElem = lenV - 1;
    constexpr std::int64_t zero64 = 0;

    if (IsRep) {
        testMax = v[lastElem] * m;
    } else if (IsMult) {
        const int lenMinusM = freqs.size() - m;

        for (int i = freqs.size() - 1; i >= lenMinusM; --i) {
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

            if (repsCounter[ind] == 0) {
                ++ind;
            }

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

    for (int i = 0; i < m; ++i) {
        finalCheck += v[z[i]];
    }

    if (finalCheck != target) {
        return 0;
    }

    return 1;
}

void GetTarget(const std::vector<double> &v,
               const std::vector<int> &Reps,
               PartDesign &part, int m, int lenV) {

    part.width = m;
    std::vector<int> z(m, 0);
    std::vector<int> zExpanded;

    std::vector<std::int64_t> v64(v.cbegin(), v.cend());
    const int lastCol = m - 1;

    for (std::size_t i = 0; i < Reps.size(); ++i) {
        for (int j = 0; j < Reps[i]; ++j) {
            zExpanded.push_back(i);
        }
    }

    const int res = GetFirstPartition(v64, zExpanded, z, Reps,
                                      part.target, m, lastCol, lenV,
                                      part.isRep, part.isMult);

    if (res == 1) {
        part.startZ    = z;
        part.solnExist = true;
        part.mapTar    = std::accumulate(z.cbegin(), z.cend(), 0) +
                         static_cast<int>(!part.includeZero) * part.width;

        if ((part.mapTar * part.slope - part.target) % (part.width)) {
            cpp11::stop("Strange mapping!!!");
        }

        part.shift = (part.mapTar * part.slope - part.target) / part.width;
    } else {
        part.startZ.assign(part.width, 0);
        part.solnExist = false;
    }
}

void SetStartPartitionZ(const std::vector<int> &Reps,
                        PartDesign &part) {

    part.startZ.assign(part.width, 0);

    switch (part.ptype) {
        case PartitionType::LengthOne: {
            part.startZ.back() = part.target;
            break;
        } case PartitionType::RepStdAll: {
            part.startZ.back() = part.target;
            break;
        } case PartitionType::RepNoZero: {
            if (part.isWeak && part.includeZero) {
                part.startZ.back() = part.target;
            } else {
                std::fill(part.startZ.begin(), part.startZ.end(), 1);
                part.startZ.back() = part.target - part.width + 1;
            }
            break;
        } case PartitionType::RepShort: {
            part.startZ.back() = part.target;
            break;
        } case PartitionType::DstctStdAll: {
            part.startZ.back() = part.target;
            break;
        } case PartitionType::DstctNoZero: {
            std::iota(part.startZ.begin(), part.startZ.end(), 1);
            part.startZ.back() = part.target - (part.width *
                                            (part.width - 1)) / 2;
            break;
        } case PartitionType::DstctOneZero: {
            std::iota(part.startZ.begin(), part.startZ.end(), 0);
            part.startZ.back() = part.target - ((part.width - 1) *
                                            (part.width - 2)) / 2;
            break;
        } case PartitionType::DstctMultiZero: {
            if (Reps.front() >= (part.width - 1)) {
                part.startZ.back() = part.target;
            } else {
                std::iota(part.startZ.begin() + Reps.front(),
                          part.startZ.end(), 1);
                part.startZ.back() = part.target - (part.width -
                    Reps.front()) * (part.width - (Reps.front() + 1)) / 2;
            }
            break;
        } default: {
            part.startZ.back() = part.target;
        }
    }
}

int DiscoverPType(const std::vector<int> &Reps,
                  PartDesign &part, int lenV) {

    if (part.ptype == PartitionType::RepCapped) {
        std::vector<int> isoz(part.width, 0);
        isoz.back() = part.mapTar -
            static_cast<int>(!part.includeZero) * part.width;

        if (part.isWeak && isoz == part.startZ) {
            part.ptype = PartitionType::RepNoZero;
            return 1;
        } else if (part.isComp && isoz == part.startZ) {
            part.ptype = part.includeZero ? PartitionType::RepShort :
                                            PartitionType::RepNoZero;
            return 1;
        } else if (isoz == part.startZ) {
            part.ptype = PartitionType::RepNoZero;
            return 1;
        }
    } else {
        for (auto ptype: DistPTypeArr) {
            std::vector<int> isoz(part.width, 0);

            switch (ptype) {
                case PartitionType::DstctNoZero: {
                    std::iota(isoz.begin(), isoz.end(), 0);
                    isoz.back() = part.mapTar - 1 - (part.width *
                        (part.width - 1)) / 2;

                    break;
                } case PartitionType::DstctMultiZero: {
                    if (!Reps.empty() && Reps.front() >= (part.width - 1)) {
                        isoz.back() = part.mapTar;
                    } else if (!Reps.empty()) {
                        std::iota(isoz.begin() + Reps.front(),
                                  isoz.end(), 1);
                        isoz.back() = part.mapTar - (part.width -
                            Reps.front()) * (part.width - (Reps.front() + 1)) / 2;
                    }

                    break;
                } case PartitionType::DstctCappedMZ: {
                    if (!Reps.empty()) {
                        int testSum = 0;
                        int max_rep = part.width;

                        for (int j = part.cap; max_rep > 1; --j) {
                            testSum += j;
                            --max_rep;

                            if (testSum >= part.mapTar) {
                                break;
                            }
                        }

                        max_rep = std::min(max_rep, Reps.front());

                        if (max_rep < (part.width - 1)) {
                            std::iota(isoz.begin() + max_rep,
                                      isoz.end(), 1);
                        }

                        int bound = part.mapTar - (part.width - max_rep) *
                                    (part.width - (max_rep + 1)) / 2;
                        int idx = part.width - 1;

                        for (int my_max = part.cap; my_max < bound &&
                             isoz[idx]; --idx, --my_max) {

                            bound += (isoz[idx - 1] - my_max);
                            isoz[idx] = my_max;
                        }

                        isoz[idx] = bound;
                    }

                    break;
                } default: {
                    // Do nothing
                }
            }

            if (isoz == part.startZ && ptype != PartitionType::DstctCapped) {
                if (part.isMult && part.allOne) {
                    part.ptype = ptype;
                    return 1;
                } else if (!part.isMult && !part.isRep) {
                    part.ptype = ptype;
                    return 1;
                }
            }
        }
    }

    return 0;
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
// If includeZero = TRUE:
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
// Since we cannot guarantee that a will evenly divide T or the smallest
// element of v, we instead use the following form:
//
//                  T_mapped = (T - min(v) * m) / a
//
// This is true regardless of IsRep, IsMult, etc. as our original
// assumption is that we are given the ith partition.
// ****************************************************************

// Our philosophy is that if the user supplies the width, we
// give that to them. If the user leaves it NULL, we figure
// out the 'best' width possible. You will see that we
// give preference to the longest possible width.
void StandardDesign(const std::vector<int> &Reps,
                    PartDesign &part, int m, int lenV) {

    // The excerpt below only holds when elements are
    // distinct (except zero)
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

    const double discriminant = 1.0 + 8.0 * static_cast<double>(part.target);
    const int max_width       = (-1 + std::sqrt(discriminant)) / 2;

    // We set an estimated initial width if it isn't given by the user
    int width = part.isRep ? m : (part.mIsNull ? std::max(max_width, 1) : m);
    part.solnExist = true;

    if (width == 1) {
        part.ptype = PartitionType::LengthOne;
    } else if (part.isMult) {

        // We need lenV == target + 1 because we could have a case where
        // we have zero and allOne, but we are missing at least one
        // element in order to guarantee we generate all possible
        // partitions: E.g. v = 0, 1, 2, 3, 4; Reps = c(4, rep(1, 4));
        // target = 6. This is equivalent to checking v.back() == target
        if (part.includeZero && lenV == part.target + 1 && part.allOne) {

            // We need at least 1 non-zero value in
            // order to set getAll = true. E.g
            //
            // Case: Distinct partitions of 10
            // width = 4 (determined above)
            //
            // 1st part:   0 0 0 10
            //
            // Thus zero must have at least a
            // frequency of 3 = 4 - 1

            if (width == max_width) {
                if ((Reps.front() >= (max_width - 1))) {
                    part.ptype = PartitionType::DstctStdAll;
                } else {
                    part.ptype = PartitionType::DstctMultiZero;
                }
            } else if (width <= (max_width + Reps.front())) {
                part.ptype = PartitionType::DstctMultiZero;
            } else {
                part.solnExist = false;
            }
        } else {
            // N.B. In the calling function we have ensured that if the
            // freqs arg is invoked with all ones, we set IsMult to false.
            // From this point, we have standard multiset partition case
            part.ptype = PartitionType::Multiset;

            if (part.mIsNull) {
                width = 0;

                // If includeZero = true, we start at index 1
                std::int64_t testTar = 0;
                bool vBigEnough = false;

                for (int i = part.includeZero, val = 1;
                     i < lenV && testTar < part.target; ++i, ++val) {

                    testTar += Reps[i] * val;
                    width   += Reps[i];

                    if (testTar > part.target) {
                        const int quotient = (testTar - part.target
                                                  + val + 1) / val;
                        width -= quotient;
                        vBigEnough = true;
                        break;
                    }
                }

                part.solnExist = vBigEnough;
            }
        }
    } else if (part.isRep) {
        if (part.isWeak && part.includeZero) {
            if (part.mIsNull) width = part.target;
            part.ptype      = PartitionType::RepNoZero;
            part.mapTar    += width;
            part.mapIncZero = false;
        } else if (part.mIsNull && part.includeZero) {
            width      = part.target; // i.e. 1 * target = target
            part.ptype = PartitionType::RepStdAll;
        } else if (part.mIsNull) {
            width      = part.target; // i.e. 1 * target = target
            part.ptype = PartitionType::RepNoZero;
        } else if (part.isComp && part.includeZero && width < part.target) {
            part.ptype = PartitionType::RepShort;
        } else if (part.includeZero && width < part.target) {
            part.ptype      = PartitionType::RepShort;
            part.mapTar    += width;
            part.mapIncZero = false;
        } else if (!part.isComp && !part.isComb && part.includeZero) {
            part.ptype = PartitionType::RepStdAll;
        }else if (part.includeZero) {
            width      = part.target;
            part.ptype = PartitionType::RepStdAll;
        } else if (width <= part.target) {
            part.ptype = PartitionType::RepNoZero;
        } else {
            part.ptype = PartitionType::RepNoZero;
            part.solnExist = false;
        }
    } else {
        if (part.includeZero && part.isComp) {
            part.ptype = PartitionType::DstctOneZero;
        } else if (part.includeZero) {
            part.ptype = PartitionType::DstctOneZero;
             // We need to add m in target in order to
             // correctly count the number of partitions
            part.mapTar += width;
            part.mapIncZero = false;
        } else {
            part.ptype = PartitionType::DstctNoZero;
        }
    }

    part.width = width;

    // Since this is only called in the standard case, this calcuation
    // is unnecessary. That is, part.mapTar = part.target and
    // part.slope = 1 in this case, which means the numerator is
    // equivalent to zero ==>> part.shift = 0. We leave it here
    // for completeness
    part.shift = (part.mapTar * part.slope - part.target) / part.width;
}

void CheckPartition(const std::vector<std::string> &compFunVec,
                    const std::vector<double> &v, const std::string &mainFun,
                    const std::vector<double> &target, PartDesign &part,
                    int lenV, int m, double tolerance, bool IsBetween) {

    /// We start by assuming we don't have a nice partition case
    part.ptype = PartitionType::NotPartition;
    bool IsPartition = false;

    if (compFunVec.front() == "==" && mainFun == "sum") {
        if (static_cast<std::int64_t>(v[0]) == v[0]) {

            IsPartition = true;
            const double tarDiff = (v.size() > 1) ? v[1] - v[0] : 1;

            for (std::size_t i = 1; i < v.size(); ++i) {
                const double testDiff = v[i] - v[i - 1];

                // We must multiply by m (the length of our partitions) as
                // we will eventually be adding these values m times.
                if (m * std::abs(testDiff - tarDiff) > tolerance ||
                    static_cast<std::int64_t>(v[i]) != v[i]) {
                    IsPartition = false;
                    break;
                }
            }
        }

        if (IsPartition &&
            (target.size() == 1 || target.front() == target.back()) &&
            static_cast<std::int64_t>(target.front()) == target.front()) {
            part.target = target.front();
        } else {
            IsPartition = false;
        }
    }

    if (!IsPartition &&
        (compFunVec.front() == "==" || IsBetween) &&
        mainFun != "max" &&
        mainFun != "min" &&
        !part.mIsNull) {

        // N.B. When we arrive here, the user must provide the width.
        // That is, m cannot be NULL
        part.ptype = PartitionType::CoarseGrained;
    }

    part.isPart = IsPartition;
}

// Right now, we have no fast method for calculating the number of partitions
// of multisets, so the variable bIsCount, is used only when we call
// partitionCount from R. If we are actually generating results, this will
// be set to false.
void SetPartitionDesign(const std::vector<int> &Reps,
                        const std::vector<double> &v,
                        PartDesign &part, ConstraintType &ctype,
                        int lenV, int &m, bool bIsCount) {

    // Now that we know we have partitions, we need to determine
    // if we are in a mapping case. There are a few of ways
    // this can happen.
    //
    // 1. If the first element isn't zero or one.
    // 2. If the distance between elements is greater than 1.
    // 3. If the first element is zero or one and the distance
    //    between elements is 1 and the target isn't equal to
    //    the largest value in v.
    //
    // The comment starting with ***** and ending with *****
    // are here for future versions that will determine m (i.e.
    // mIsNull = TRUE) when the user doesn't supply it. For right
    // now, since we require that the elements be spaced evenly
    // (i.e cases such as 0, 3:14 will not be considered), these
    // extra checks are superfulous.
    //
    // *****
    // For point 2, we need to check the difference between two
    // elements of v. If the first element is 0, we have satisfied
    // the first point, and if the next 10 elements are
    // seq(1, 99, 10), then only checking the difference between
    // the first 2 elements would not reveal the nonstandard
    // distance between each element.
    // *****
    //
    // Note, we have already ensured above that if we have
    // negative values and the differenece between every value
    // isn't the same, then we don't meet the partition scenario.

    part.slope = (v.size() > 1) ? v[1] - v[0] : 1;

    // When allOne = true with isMult = true, we can apply optimized
    // algorithms for counting distinct partitions of varying lengths.
    part.allOne = part.isMult ? std::all_of(Reps.cbegin() + 1,
                                            Reps.cend(), [](int v_i) {
                                                return v_i == 1;
                                            }) : false;

    // When we have repetition or the distinct case, part.isMult will
    // be false as well as part.allOne will be false. When part.isMult is
    // true, the only way we are in the standard case is when part.allOne
    // is also true. Thus the expression below.
    const bool standard_freq = part.isMult == part.allOne;
    const bool zero_or_one   = v.front() == 0 || v.front() == 1;
    const bool standard_dist = part.slope == 1;
    const bool standard_tar  = v.back() == part.target;

    if (zero_or_one && standard_dist && standard_tar && standard_freq) {
        // Remember, lenV is the length of the vector v, so we could have a
        // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
        // cause a problem if we were to allow this. We have already ensured
        // that the distance between each element is the same. This means for
        // the example we gave, we would have length(unique(diff(v))) > 1,
        // which means PartitionCase would be false, and thus the general
        // algorithm would be executed.

        part.includeZero = (v.front() == 0);
        part.mapIncZero  = part.includeZero;
        part.isWeak      = part.isWeak && part.includeZero;
        part.mapTar      = part.target;
        part.cap         = v.back();

        ctype = ConstraintType::PartStandard;
        StandardDesign(Reps, part, m, lenV);
        SetStartPartitionZ(Reps, part);
    } else {
        // For right now, if m is not provided (i.e. mIsNull = true),
        // we don't try to figure out the appropriate length. Note,
        // this only applies to non-canonical partitions.
        part.mIsNull = false;

        // We can only have weak compositions when zero is included. We
        // can't use part.includeZero for the reasons in the below comment
        part.isWeak = part.isWeak && (v.front() == 0);

        // When we are mapping cases, it is easy to calculate the number of
        // results when we map zero to one, since the mapped value  will
        // count zero as one when weak = TRUE.
        part.includeZero = part.allOne || (part.isComp && v.front() == 0 && !part.isWeak);
        part.mapIncZero  = part.includeZero;
        part.cap         = lenV - part.mapIncZero;

        part.ptype = (m == 1) ? PartitionType::LengthOne :
            (part.isMult ? PartitionType::Multiset :
            (part.isRep ? PartitionType::RepCapped :
                 PartitionType::DstctCapped));

        ctype = ConstraintType::PartMapping;
        GetTarget(v, Reps, part, m, lenV);

        if (part.solnExist && part.ptype != PartitionType::LengthOne) {
            DiscoverPType(Reps, part, lenV);
        }
    }

    PartitionsCount(Reps, part, lenV, bIsCount);
}
