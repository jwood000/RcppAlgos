
// void BinaryNextElem(int &uppBnd, int &lowBnd, int &ind, int lastElem,
//                     std::int64_t target, std::int64_t partial,
//                     const std::vector<std::int64_t> &v) {
//     
//     std::int64_t dist = target - (partial + v[ind]);
//     
//     while ((uppBnd - lowBnd) > 1 && dist != 0) {
//         const int mid = (uppBnd - lowBnd) / 2;
//         ind = lowBnd + mid;
//         dist = target - (partial + v[ind]);
//         
//         if (dist > 0) {
//             lowBnd = ind;
//         } else {
//             uppBnd = ind;
//         }
//     }
//     
//     // Check last index. N.B. There are some cases when ind == lowBnd and dist < 0.
//     // This will not matter as we simply reassign ind and recompute dist 
//     if (dist < 0) {
//         ind = lowBnd;
//         dist = target - (partial + v[ind]);
//     }
//     
//     // We must have dist <= 0. Below is an informal proof.
//     // The sub-sequences are defined as below:
//     //                  A_max = a_(i + 1), a_(i + 2), ..., a_m
//     //                  A_set = a_1, a_2, ..., a_(i - 1)
//     // A_set are those elements that have already been determined by the algorithm.
//     // A_max is maximal (i.e. constructed of the the (m - i) largest elements). We
//     // seek to determine the i_th element given the following contraints:
//     //                      A_sum = A_set + a_i + A_max
//     //                       dist = target - A_sum
//     // With the goal of finding the minimum lexicographic combination such that the
//     // dist = 0 (i.e. target = A_sum). If we have dist > 0 for any i, then it will
//     // be impossible to obtain dist = 0. dist > 0 implies that the target > A_sum,
//     // and since A_max is already maximal, we are not able to increase A_sum in
//     // later iterations, thus we must have dist <= 0 for all i.
//     
//     if (dist > 0 && ind < lastElem) {
//         ++ind;
//     }
// }
// 
// int GetFirstPartition(int m, const std::vector<std::int64_t> &v,
//                       bool IsRep, bool IsMult, std::vector<int> &z,
//                       const std::vector<int> &freqs, std::int64_t target,
//                       std::vector<int> repsCounter, int lastCol, int lastElem) {
//     
//     std::int64_t testMax = 0;
//     constexpr std::int64_t zero64 = 0;
//     
//     if (IsRep) {
//         testMax = v[lastElem] * m;
//     } else if (IsMult) {
//         const int lenMinusM = freqs.size() - m;
//         
//         for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j){
//             testMax += v[freqs[i]];
//         }
//         
//     } else {
//         testMax = std::accumulate(v.cend() - m, v.cend(), zero64);
//     }
//     
//     if (testMax < target)  {return 0;}
//     int zExpCurrPos = IsMult ? freqs.size() - m : 0;
//     int currPos = IsMult ? freqs[zExpCurrPos] : (IsRep ? lastElem : (v.size() - m));
//     
//     std::int64_t partial = testMax;
//     partial -= v[currPos];
//     std::int64_t testMin = 0;
//     
//     if (IsRep) {
//         testMin = v[0] * m;
//     } else if (IsMult) {
//         
//         for (int i = 0; i < m; ++i) {
//             testMin += v[freqs[i]];
//         }
//         
//     } else {
//         testMin = std::accumulate(v.cbegin(), v.cbegin() + m, zero64);
//     }
//     
//     if (testMin > target) {return 0;}
//     int mid = currPos / 2;
//     std::int64_t dist = target - (partial + v[mid]);
//     
//     int lowBnd = (dist > 0) ? mid : 0;
//     int uppBnd = (dist > 0) ? currPos : mid;
//     int ind = mid;
//     
//     for (int i = 0; i < lastCol; ++i) {
//         BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
//         z[i] = ind;
//         partial += v[ind];
//         
//         if (IsMult) {
//             --repsCounter[ind];
//             
//             if (repsCounter[ind] == 0)
//                 ++ind;
//             
//             ++zExpCurrPos;
//             currPos = freqs[zExpCurrPos];
//         } else if (!IsRep) {
//             ++ind;
//             ++currPos;
//         }
//         
//         lowBnd = ind;
//         uppBnd = currPos;
//         mid = (uppBnd - lowBnd) / 2;
//         
//         ind = lowBnd + mid;
//         partial -= v[currPos];
//     }
//     
//     BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
//     z[lastCol] = ind;
//     
//     // The algorithm above finds the first possible sum that equals
//     // target. If there is no combination of elements from v that sum
//     // to target, the algo returns the combination such that its sum
//     // is closest to target and greater than target
//     std::int64_t finalCheck = 0;
//     
//     for (int i = 0; i < m; ++i)
//         finalCheck += v[z[i]];
//     
//     if (finalCheck != target)
//         return 0;
//     
//     return 1;
// }
// 
// template <typename T>
// int GetMappedTarget(const std::vector<T> &v, std::int64_t target,
//                     const std::vector<int> &Reps, int m, int lenV,
//                     bool IsMult, bool IsRep) {
//     
//     std::vector<int> z(m, 0);
//     std::vector<int> zExpanded;
//     int myTarget = 0;
//     
//     std::vector<std::int64_t> v64(v.cbegin(), v.cend());
//     const int lastCol = m - 1;
//     
//     for (std::size_t i = 0; i < Reps.size(); ++i)
//         for (int j = 0; j < Reps[i]; ++j)
//             zExpanded.push_back(i);
//     
//     Rcpp::Rcout << "\n";
//     Rcpp::print(Rcpp::wrap(target));
//     Rcpp::print(Rcpp::wrap(lastCol));
//     Rcpp::print(Rcpp::wrap(lenV - 1));
//     Rcpp::print(Rcpp::wrap(m));
//     Rcpp::print(Rcpp::wrap(v64));
//     
//     if (GetFirstPartition(m, v64, IsRep, IsMult, z, zExpanded,
//                           target, Reps, lastCol, lenV - 1)) {
//         
//         Rcpp::Rcout << "\n";
//         Rcpp::print(Rcpp::wrap("First combo index"));
//         Rcpp::print(Rcpp::wrap(z));
//         
//         myTarget = std::accumulate(z.begin(), z.end(), 0) + m;
//     }
//     
//     return myTarget;
// }
// 
// template <typename T>
// void GetPartitionCase(const std::vector<std::string> &compFunVec, std::vector<T> &v,
//                       const std::string &mainFun, const std::vector<T> &target,
//                       Sign mySign, PartitionType &PartType, ConstraintType &ConstType,
//                       Sign mySign, distinctType &distinctTest, SEXP Rlow,
//                       std::vector<int> &Reps, int lenV, int &m, double tolerance,
//                       bool IsMult, bool IsRep, bool IsBetween, bool mIsNull) {
//     
//     ConstType = ConstraintType::General;
//     bool bLower = false;
//     
//     // Currently, we are not able to generate the nth
//     // lexicographical partition. Thus, if lower is
//     // non-trivial, we must use most general algo.
//     if (!Rf_isNull(Rlow)) {
//         auto tempLower = FromCpp14::make_unique<mpz_t[]>(1);
//         mpz_init(tempLower[0]);
//         
//         createMPZArray(Rlow, tempLower.get(), 1, "lower");
//         bLower = mpz_cmp_si(tempLower[0], 1) > 0;
//     }
//     
//     // compFunVec should be non-empty if we made it this far.
//     // Doesn't hurt to check
//     if (!compFunVec.empty() && !bLower) {
//         
//         /// We start by assuming we don't have a nice partition case
//         bool PartitionCase = false;
//         std::int64_t tarTest = 0;
//         
//         if (compFunVec[0] == "==" && mainFun == "sum" && lenV >  1 && m >  1) {
//             // We need to make sure not to include zero in the check below.
//             // This is so because the zero can be used to obtain partitions
//             // of differing lengths. Under normal circumcstances, this is
//             // no problem because we simply have 0, 1, 2,..., however with
//             // the capped cases (i.e. they don't start at 0 or 1, e.g. 3:14)
//             // this case would be excluded because (3 - 0) != (4 - 3). Note,
//             // We can only do this if all values have the same sign. When we
//             // have mixed signs numbers, indexing breaks down. We are making
//             // the case that we can't get all partititions of every length
//             // when this occurs. First off, we have mapping issues. E.g. Let
//             // v: -15, -9, -3, 0, 3, 9,..., 99 and a target of 93. For m = 3,
//             // the first partition is 0, 0, -15, 9, 99 which maps to 
//             // 4, 4, 1, 6, 21 giving a new target of 37. Now let m = 5
//             // The first partition is -15, -15, -15, 39, 99 which maps to 
//             // 1, 1, 1, 11, 21 for a mapped target of 35 (which is not 37!).
//             //
//             // Secondly, even if we could map properly for differing lengths
//             // we would have issues with ordering (lexicographically).
//             
//             std::vector<T> pTest;
//             
//             if (mySign != Sign::MixedBag) {
//                 for (auto val: v) {
//                     if (val != 0) {
//                         pTest.push_back(val);
//                     }
//                 }
//             } else {
//                 pTest = v;
//             }
//             
//             std::sort(pTest.begin(), pTest.end());
//             const T tarDiff = pTest[1] - pTest[0];
//             
//             if (static_cast<std::int64_t>(pTest[0]) == pTest[0]) {
//                 PartitionCase = true;
//                 
//                 for (std::size_t i = 1; i < pTest.size(); ++i) {
//                     const T testDiff = pTest[i] - pTest[i - 1];
//                     
//                     if (std::abs(testDiff - tarDiff)  > tolerance ||
//                         static_cast<std::int64_t>(pTest[i]) != pTest[i]) {
//                         
//                         PartitionCase = false;
//                         break;
//                     }
//                 }
//             }
//             
//             if (target.size() == 1 || target.front() == target.back()) {
//                 tarTest = static_cast<std::int64_t>(target.front());
//                 
//                 if (PartitionCase)
//                     PartitionCase = (tarTest == target.front());
//             } else {
//                 PartitionCase = false;
//             }
//         }
//         
//         if (PartitionCase) {
//             // Now that we know we have partitions, we need to determine
//             // if the final result needs to be mapped. There are a couple
//             // of ways this can happen.
//             // 
//             // 1. If the first element isn't zero or one.
//             // 2. If the distance between elements is greater than 1.
//             //
//             // The vector vBase will take on the underlying base partition.
//             //
//             // Note, we have already ensured above that if we have
//             // negative values and the differenece between every value
//             // isn't the same, then we don't meet the partition scenario.
//             
//             std::vector<int> vBase(v.size());
//             int mappedTarget = 0;
//             const T testDiff = v[1] - v[0];
//             const bool condition_1 = (v.front() != 1 && v.front() != 0) || testDiff != 1;
//             
//             if (condition_1 && !Sign::MixedBag) {
//                 // Set our constraint type to indicate mapping
//                 // will be needed. See PartitionMain.cpp
//                 ConstType = ConstraintType::PartMapping;
//                 std::iota(vBase.begin(), vBase.end(), 0);
//                 mappedTarget = GetMappedTarget(v, tarTest, Reps,
//                                                m, lenV, IsMult, IsRep);
//             } else {
//                 ConstType = ConstraintType::PartStandard;
//                 mappedTarget = tarTest;
//                 
//                 for (int i = 0; i < lenV; ++i)
//                     vBase[i] = v[i];
//             }
//             
//             Rcpp::Rcout << "\n";
//             Rcpp::print(Rcpp::wrap("mappedVector"));
//             Rcpp::print(Rcpp::wrap(vBase));
//             Rcpp::Rcout << "\n";
//             Rcpp::print(Rcpp::wrap("mappedTarget"));
//             Rcpp::print(Rcpp::wrap(mappedTarget));
//             Rcpp::stop("digg");
//             
//             // We sorted v above to ensure that the last element is the maximum
//             const int myMax = vBase.back();
//             const int IncludeZero = (vBase.front() == 0);
//             
//             // Remember, lenV is the length of the vector v, so we could have a
//             // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
//             // cause a problem if we were to allow this. We have already ensured 
//             // that the distance between each element is the same. This means for
//             // the example we gave, we would have length(unique(diff(v))) > 1,
//             // which means PartitionCase would be false, and thus the general
//             // algorithm would be executed.
//             //
//             // We do have to ensure that the smallest element is non-negative, othe-
//             // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & limit = 10
//             // would pass as v = 0:9, m = 7, rep = TRUE, & limit = 9, --or--
//             // v = 1:10, m = 7, rep = TRUE, & limit = 10 (Hence v.front() >= 0)
//             
//             if (                
//                     myMax == mappedTarget
//                 && lenV + IncludeZero == mappedTarget
//                 && vBase.front() >= 0
//             )
//             {
//                 distinctTest = DistinctAttr(lenV, m, IsRep, IsMult, mappedTarget,
//                                             Reps, static_cast<bool>(IncludeZero), mIsNull);
//                 
//                 if (distinctTest.limit > 0) {
//                     m = distinctTest.limit;
//                     
//                     if (IncludeZero) {
//                         if (IsMult) {
//                             if (distinctTest.getAll) {
//                                 PartType = PartitionType::DstctStdAll;
//                             } else if (Reps[0] >= (m - 1)) {
//                                 PartType = PartitionType::DstctShort;
//                             } else {
//                                 PartType = PartitionType::DstctSpecial;
//                             }
//                         } else {
//                             PartType = PartitionType::DstctOneZero;
//                         }
//                     } else {
//                         PartType = PartitionType::DstctNoZero;
//                     }
//                 } else if (IsRep) {
//                     if (IncludeZero) {
//                         PartType = PartitionType::Traditional;
//                     } else {
//                         PartType = PartitionType::TradNoZero;
//                     }
//                     
//                     if (m >= lenV) { 
//                         if (IncludeZero) {
//                             m = lenV - 1;
//                         } else {
//                             ConstType = ConstraintType::General;
//                         }
//                     }
//                 }
//             }  else {
//                 if (IsRep) {
//                     PartType = PartitionType::TradCapped;
//                 } else if (!IsMult) {
//                     PartType = PartitionType::DistCapped;
//                 }
//             }
//         } else if (
//                 (compFunVec[0] == "==" || IsBetween)
//             && lenV > 1
//         && m > 1
//         && mainFun != "max"
//         && mainFun != "min"
//         )
//         {
//             
//             // N.B. When we arrive here, the user must provide the width.
//             // That is, m cannot be NULL
//             ConstType = ConstraintType::PartitionEsque;
//         }
//     }
// }