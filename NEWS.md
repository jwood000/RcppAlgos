# RcppAlgos 2.8.2

## New Features:

* We have added `comboGroupsIter`. This a flexible iterator for iterating over partitions of groups. Offers random access, the ability to retrieve the next group or the next `n` groups at a time while keeping memory low. 

# RcppAlgos 2.8.1

## Bug Fixes:

* Fixed integer overflow bug when converting vector size to number of rows in a matrix. See [Issue #45](<https://github.com/jwood000/RcppAlgos/issues/45>).

* Fixed constraint permutation iterator. Before this fix, if a user requested a certain number of results such that the `std::next_permutation` algorithm was not exhausted, the next iteration requested would be nonsensical.

# RcppAlgos 2.8.0

## New Features:

* `comboGroups` can now handle groups of different sizes.

# RcppAlgos 2.7.2

## Other:

* Using a copy of `gmpxx.h` source in order to easily build on all platforms.

# RcppAlgos 2.7.1

## Other:

* Updated default C++ specification.

# RcppAlgos 2.7.0

## Other:

* Now using `gmpxx.h` instead of `gmp.h`.

# RcppAlgos 2.6.0

## New Features:

* Added integer composition functions: `compositionsCount`, `compositionsGeneral`, `compositionsSample`, `compositionsIter`, and `compositionsRank`.

# RcppAlgos 2.5.5

## New Features:

* Added ranking functions: `comboRank`, `permuteRank`, and `partitionsRank`.

* Added back the ability to interrupt general constraint problems via `cpp11::check_user_interrupt()`.

## Other:

* Corrected `if()` conditions comparing `class()` to string. Changed to `is.character`.

## Bug Fixes:

* Now checking class of input vector for partition funcitons.

* Now when `partitionsCount` returns 0, the number of results is zero. Before, we were checking for count of partitions to be greater than zero, otherwise we would use the standard combinatorial counting functions to determine the number of results. This lead to strange results with elements not present in the original vector.

* For `partitionsSample`, in cases when we would rely on generating the partitions one at a time until there are no more (e.g. with `partitionsGeneral`), the number of partitions isn't calculated. This leads to the error: "n exceeds the maximum number of possible results". This is now fixed.

# RcppAlgos 2.5.4

## Other:

* Added missing includes

# RcppAlgos 2.5.3

## Other:

* Fixed urls

# RcppAlgos 2.5.2

## Bug Fixes:

* Fixed valgrind memory issue with package `cpp11`. Now using `cpp11::stop` instead of `Rf_error` so as to avoid `longjmp` (the root cause of the memory issues).

# RcppAlgos 2.5.1

## Bug Fixes:

* Fixed memory issue when the number of results under contraint is less than the requested number.

# RcppAlgos 2.5.0

## New Features:

* Added partition specific functions: `partionsGeneral`, `partitionsCount`, `partitionsIter`, and `partitionsSample`

* `comboIter` and `permuteIter` now work with constraints.

* Dropped `Rcpp` and `RcppThread` as a dependency to reduce compile time and binary size. As a result, there is no longer the ability to interrupt long running processes. Will investigate in next release.

* We have added the parameter `FUN.VALUE` to `comboGeneral` and `permuteGeneral`. It acts as a template for the return value of `FUN`. The implementation is modeled after `vapply`.

## Enhancements:

* Improved underlying algorithm for `comboGrid` to be more memory efficient.

* Made minor changes to make data types more consistent in `primeCount` and `primeSieve`.

## Bug Fixes:

* When `permuteGeneral` is used with multisets and the width is maximized, multithreading would fail. This is fixed in 2.5.0.

* Fixed bug in retreiving the nth result in `comboGroup` and `comboGroupSample` when the number of results was greater than `2^31 - 1` and less than `2^53 - 1`. E.g. `comboGroupsSample(27, 9, seed = 4, sampleVec = 1606990240475839)` gives incorrect results in the 5th group in prior versions. Now fixed!.

## Other:

* In this version we no longer output lexicographical composititions in very special circumstances outlined in older documentation using `permuteGeneral`. This was done for consistency as we felt that the output diverged too much from the general constrained output of `permuteGeneral` (See [Output Order with permuteGeneral](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#output-order-with-permutegeneral>)). Research is in the initial stage that is focused on implementing a new family of functions similar to the partition family of functions, but for compositions.

# RcppAlgos 2.4.3

## Other:

* Fixed old urls

# RcppAlgos 2.4.2

## New Features:

* Added function `comboGrid` for efficiently generating the Cartesian product where order does not matter.

## Enhancements:

* Refactored code base to reduce binary size

## Other:

* Removed LazyData from DESCRIPTION. Also added `rmarkdown` to Suggests.

# RcppAlgos 2.4.1

## Bug Fixes:

* Removed broken URL in NEWS.Rd

# RcppAlgos 2.4.0

## New Features:

* Added `comboIter` and `permuteIter`. These functions return iterators for iterating over combinations and permutations. They have a similar interface to `comboGeneral` and `permuteGeneral` and currently only work with standard combinations and permutations. They do not yet work with constraints (This will be the focus of the next release).

* Added "High Performance Benchmarks" and "Combinatorial Iterators in RcppAlgos" vignettes

## Enhancements:

* Now able to interrupt general constraint problems (See [Interrupt Execution with Rcpp::checkUserInterrupt](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#safely-interrupt-execution-with-rcppcheckuserinterrupt>))

## Bug Fixes:

* In 2.3.5 and 2.3.6, we mistakenly allowed a `constraintFun` to be applied to a logical vector which was crashing R. We have corrected this in 2.4.0.

* Changed the data type for sizing the index matrix in `permuteGeneral`. Originally, we were using `int` and when the output was large enough, it was causing an integer overflow thus causing the index matrix to be improperly sized. We have since changed the data type to the recommended `std::size_t` (See [Is there a max array length limit in C++?](<https://stackoverflow.com/q/216259/4408538>))

# RcppAlgos 2.3.6

## Bug Fixes:

* Fixed bug associated with integer vectors, multisets, and constraints. See [Issue #12](<https://github.com/jwood000/RcppAlgos/issues/12#issue-467908452>) for more information.

# RcppAlgos 2.3.5

## New Features:

* Added `comboGroups`, `comboGroupsCount`, and `comboGroupsSample`. These functions deal with partitioning a vector/set into groups of equal size. See [Combinations in R by Groups](<https://stackoverflow.com/q/57732672/4408538Create>). See the related integer sequences A025035-A025042 at https://oeis.org (E.g. https://oeis.org/A025036A025036 for Number of partitions of `(1, 2, ..., 4n)` into sets of size 4.)

* Added vignettes (First version with vignettes)

* Added website via the excellent package `pkgdown`

* Now using C++11 instead of C++14. See [Issue #10](<https://github.com/jwood000/RcppAlgos/issues/10>) for more information.

## Enhancements:

* Extended general partitions algorithm to multisets. E.g. `comboGeneral(10, 8, freqs = rep(1:5, 2), constraintFun = "sum", comparisonFun = "==", limitConstraints = 55)`

* Improved constraint algorithm for the general case.

* Added support for complex and raw types for all combinatorial functions.

* Improved permutation algorithm for all cases (I.e. no rep, with rep, and with multisets).

* Added loop unrolling to prime sieve algorithm for improved efficiency.

## Bug Fixes:

* Corrected checks for total number of partitions and assignment of number of rows when `upper` is applied in `{combo|permute}General`. See [Issue #9](<https://github.com/jwood000/RcppAlgos/issues/9#issue-467908452Issue%20#9>) for more information.

* `permuteGeneral` no longer alters source vector. See [Issue #11](<https://github.com/jwood000/RcppAlgos/issues/11>) for more information.

# RcppAlgos 2.3.4

* Fixed clang/gcc-ASAN and valgrind issues in `2.3.3`. These issues were arising from finding the first vector to meet the criteria in PartitionRep/Distinct. We also found further issues in the standard functions when the length of the partition was 2. The algorithm would eventually try and access an element of a vector at index -1. These fixes were confirmed by successful Rdevel CMD check under 'r-devel-ubsan-clang' via docker using the advice found here: https://knausb.github.io/2017/06/reproducing-a-clang-ubsan-issue/ and http://dirk.eddelbuettel.com/code/sanitizers.html.

# RcppAlgos 2.3.3

* Fixed clang-UBSAN issue in `2.3.2`. It was caused by populating a vector of ints with values larger than `2^31 - 1`.

* Added optimized algorithm to `{combo|permute}General` when `constraintFun = "sum"`, `comparisonFun = "=="`, and the vector passed has a special form. This problem is a special case of the subset sum problem.

* Using `std::vector` and `push_back` member function instead of pre-allocating matrix when constraint is applied in `{combo|permute}General`. This alleviates the need to guess the upper limit and subsequently subset as only elements that meet the constraints are added.

* Fixed error in `PollardRho.cpp` when number passed had factors close to the limit in the predefined lookup table (i.e. `constexpr int64_t FirstOmittedPrime = 3989`)

# RcppAlgos 2.3.2

* Fixed clang-UBSAN issue in `2.3.1`. It was caused by casting extremely large values to `int64_t`.

* Corrected handling of small values in `PrimeSieveBig`.

# RcppAlgos 2.3.1

* Explicitly casted to double for sqrt to silence Solaris

* Corrected handling of `NaNs`

# RcppAlgos 2.3.0

* All functions now have parallel capabilities via `RcppThread`.

* Utilizes `RMatrix.h` from `RcppParallel` for thread safe matrix class.

* Major overhaul of primeSieve for large primes.

* Added `stdThreadMax` for obtaining the number of threads available on a machine

# RcppAlgos 2.2.0

* Disabled `Parallel` argument as it was causing unpredictable errors on certain platforms. Research is ongoing to correct this for use in future versions. The development version will retain this feature.

* Corrected `UBSAN` error that caused by filling a vector of integers with numbers larger than `2^31 - 1`.

# RcppAlgos 2.1.0

* Added argument `Parallel` to general and sampling functions for increased gains in efficiency.

    * `comboGeneral(30, 10, Parallel = TRUE)`

    * `permuteGeneral(12, 7, TRUE, constraintFun = "sum", Parallel = TRUE)`

* Logical class is now preserved in combinatorial functions

* Added gmp support to combinatorial functions. Now, one can accurately and quickly work with combinations/permutations of large vectors _E.g._:

    * `comboSample(runif(10000), 100, n = 10, seed = 42, Parallel = TRUE)`

    * `permuteGeneral(factor(state.name), 20, lower = 1e15, upper = 1e15 + 1000)`

* Added `FUN` argument to all combinatorial functions. Allows user to pass custom functions to be applied to combinations/permutations.

# RcppAlgos 2.0.3

* Corrected clang `UBSAN` error identified by two different unit tests. In both situations, the problem was occurring as a result of populating a vector of integers with values from a vector of doubles that contained a NaN (Not-a-Number). Most information was obtained from Brian J. Knaus's blog titled : "Reproducing a clang-UBSAN issue" (<https://knausb.github.io/2017/06/reproducing-a-clang-ubsan-issue/>)

# RcppAlgos 2.0.2

* Corrected divide by zero in `divisorsRcpp` unit test.

# RcppAlgos 2.0.1

* Corrected spelling in `DESCRIPTION`

# RcppAlgos 2.0.0

* Changed max value and explicitly casted a few values to `int64_t` in PollardRho.cpp for efficiency while still maintaining accuracy from `2^60` to `2^62` (`isPrimeRcpp` is roughly 10% faster now).

* Updated core permutation algorithm for greater efficiency and generality

* Added capability of generating specific combinations/permutations

* Changed arguments to `comboGeneral/permuteGeneral`. `rowCap` is now `upper`.

* Added `comboSample`, `permuteSample`, `comboCount`, and `permuteCount`

* Fixed bug with `numDivisorSieve` and `divisorSieve` when the lower bound was greater than 1 and less than the `sqrt` of the upper bound. In the previous version, the numbers in this range would have duplicated values/counts.

* Increased efficiency of `numDivisorSieve` by a factor of 2.

* Updated unit tests for greater coverage. See the function `package_coverage` from the package `covr`.

# RcppAlgos 1.0.1

* Corrected precision limits in documentation from `2^64` to `2^63`.

* Changed const type in `PollardRho.cpp` from `int64_t` to `double` to correct "UndefinedBehaviorSanitizer"

* Changed examples in `primeFactorizeSieve` to reduce check time

* Added `RcppAlgos-package` man file.

# RcppAlgos 1.0.0

* Added the following functions: `primeFactorize` (vectorized pollard rho factorization), `divisorsRcpp` (vectorized factoring (complete)), `isPrimeRcpp` (vectorized primality testing using Miller-Rabin algorithm), & `primeCount` (based on the primecount algorithm by Kim Walisch)

* Completely revamped the `primeSieve` function. It is now a segmented sieve of Eratosthenes with wheel factorization based on primesieve by Kim Walisch.

* Renamed `divisorsList` to `divisorsSieve` (reason for the major version update to `1.0.0`)

* Renamed `primeFactorizationList` to `primeFactorizeSieve`

* Made the sieving functions more flexible. They are now able to generate results over a range and can also produce named objects.

* All number theoretic functions have been made more efficient. Some make use of the fast integer division library `libdivide` by ridiculousfish.

# RcppAlgos 0.2.5

* Added unit tests.

* Removed unnecessary files.

* Fixed bug in primeSieve that occurred when a number with a decimal was passed (e.g. 2.01).

* Adjusted accepted lower bound for `numDivisorSieve`, `eulerPhiSieve`, `divisorsList`, `primeSieve`, and `primeFactorizationList`.

* Fixed bug when non-unique elements are present with factors.

# RcppAlgos 0.2.4

* Fixed bug that occurs when non-unique elements are present for combinations with replacement.

# RcppAlgos 0.2.3

* Fixed segmentation fault error highlighted by valgrind check in version `0.2.2`.

* Updated `DESCRIPTION` file.

# RcppAlgos 0.2.2

* Fixed bug in constraint functions that occurred when `m = 1` and the constraint limit was equal to the last element in `v`. It was returning a 2x1 matrix with the same value twice.  It is now correctly returning a 1x1 matrix with the correct value 1 time.

* Reorganized source code such that all utility functions for the combinatoric functions are now in their own file. Additionally added header for this file.

* All combinatoric functions can now utilize the rowCap argument. Before, rowCap only applied to the combinatorial functions with constraints.

* `comboGeneral` can now find all combinations of multisets.

* Both `comboGeneral` and `permuteGeneral` can utilize the argument `m` when dealing with multisets. Before, `permuteGeneral` would simply return all permutations of a multiset. Now you can specify the lengths of the output.

# RcppAlgos 0.2.1

* Fixed bug that would occur in two edge cases involving the constraint functions.

* Slightly modified formatting for `primeSieve.Rd`

# RcppAlgos 0.2.0

* Updated combination algorithms. They are now more than twice as fast.

* Updated constraint functions so that memory access is always within container bounds

* Consolidated redundant code for outputting different `Rcpp` types (e.g. `IntegerMatrix`, `CharacterMatrix`, etc.) via a templated approach.

* Added the function permuteGeneral that is analogous to comboGeneral only instead of combinations, it gives all permutations. It has an additional argument (i.e. 'freqs') that is used to generate permutations of multisets.

* All combinatoric functions now support factor types.

# RcppAlgos 0.1.2

* Corrected minor typo in `README` file.

* Fixed minor error regarding explicitly comparing variables to large numbers that are typed out. Simply adding a decimal along with a zero remedies the situation.

# RcppAlgos 0.1.1

* Improved `ComboConstraint` function by removing unnecessary subsetting.

* Improved `PrimeSieve` internal `C++` algorithm.

* Corrected the errors with respect to the math functions in `C++`. Explicitly overloaded the parameters of these functions by casting them to the `double` type.

# RcppAlgos 0.1.0

* Initial Release
