\name{NEWS}
\title{News for Package \pkg{RcppAlgos}}

\section{Changes in RcppAlgos version 2.7.2 (Release date: 2023-02-11)}{
    \itemize{
        \item Others:
        \itemize{
            \item Using a copy of \code{gmpxx.h} source in order to easily build on all platforms.
        }
    }
}

\section{Changes in RcppAlgos version 2.7.1 (Release date: 2023-02-06)}{
    \itemize{
        \item Others:
        \itemize{
            \item Updated default C++ specification.
        }
    }
}

\section{Changes in RcppAlgos version 2.7.0 (Release date: Never Released)}{
    \itemize{
        \item Others:
        \itemize{
            \item Now using \code{gmpxx.h} instead of \code{gmp.h}.
        }
    }
}

\section{Changes in RcppAlgos version 2.6.0 (Release date: 2022-08-15)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added integer composition functions: \code{compositionsCount}, \code{compositionsGeneral}, \code{compositionsSample}, \code{compositionsIter}, and \code{compositionsRank}.
        }
    }
}

\section{Changes in RcppAlgos version 2.5.5 (Release date: Never Released)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added ranking functions: \code{comboRank}, \code{permuteRank}, and \code{partitionsRank}
            \item Added back the ability to interrupt general constraint problems via \code{cpp11::check_user_interrupt()}.
        }
        \item Other:
        \itemize{
            \item Corrected \code{if()} conditions comparing \code{class()} to string. Changed to \code{is.character}.
        }
        \item Bug Fixes:
        \itemize{
            \item Now checking class of input vector for partition funcitons.
            \item Now when \code{partitionsCount} returns 0, the number of results is zero. Before, we were checking for the count of partitions to be greater than zero, otherwise we would use the standard combinatorial counting functions to determine the number of results. This lead to strange results with elements not present in the original vector.
            \item For \code{partitionsSample}, in cases when we would rely on generating the partitions one at a time until there are no more (e.g. with \code{partitionsGeneral}), the number of partitions isn't calculated. This leads to the error: "n exceeds the maximum number of possible results". This is now fixed.
        }
    }
}

\section{Changes in RcppAlgos version 2.5.4 (Release date: Never Released)}{
    \itemize{
        \item Other:
        \itemize{
            \item Added missing includes
        }
    }
}

\section{Changes in RcppAlgos version 2.5.3 (Release date: 2022-03-31)}{
    \itemize{
        \item Other:
        \itemize{
            \item Fixed urls
        }
    }
}

\section{Changes in RcppAlgos version 2.5.2 (Release date: Never Released)}{
    \itemize{
        \item Bug Fixes:
        \itemize{
            \item Fixed valgrind memory issue with package cpp11. Now using \code{cpp11::stop} instead of \code{Rf_error} so as to avoid \code{longjmp} (the root cause of the memory issues).
        }
    }
}

\section{Changes in RcppAlgos version 2.5.1 (Release date: Never Released)}{
    \itemize{
        \item Bug Fixes:
        \itemize{
            \item Fixed memory issue when the number of results under contraint is less than the requested number.
        }
    }
}

\section{Changes in RcppAlgos version 2.5.0 (Release date: 2022-03-16)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added partition specific functions: \code{partionsGeneral}, \code{partitionsCount}, \code{partitionsIter}, and \code{partitionsSample}
            \item \code{comboIter} and \code{permuteIter} now work with constraints.
            \item Dropped Rcpp and RcppThread as a dependency to reduce compile time and binary size. As a result, there is no longer the ability to interrupt long running processes. Will investigate in next release.
            \item We have added the parameter \code{FUN.VALUE} to \code{comboGeneral} and \code{permuteGeneral}. It acts as a template for the return value of \code{FUN}. The implementation is modeled after \code{vapply}.
        }
        \item Enhancements:
        \itemize{
            \item Improved underlying algorithm for \code{comboGrid} to be more memory efficient.
            \item Made minor changes to make data types more consistent in \code{primeCount} and \code{primeSieve}.
        }
        \item Bug Fixes:
        \itemize{
            \item When \code{permuteGeneral} is used with multisets and the width is maximized, multithreading would fail. This is fixed in 2.5.0.
            \item Fixed bug in retreiving the \eqn{n^{th}}{nth} result in \code{comboGroup} and \code{comboGroupSample} when the number of results was greater than \eqn{2^{31} - 1}{2^31 - 1} and less than \eqn{2^{53} - 1}{2^53 - 1}. \emph{E.g.} \code{comboGroupsSample(27, 9, seed = 4, sampleVec = 1606990240475839)} gives incorrect results in the \eqn{5^{th}}{5th} group in prior versions. Now fixed!.
        }
        \item Other:
        \itemize{
            \item In this version we no longer output lexicographical composititions in very special circumstances outlined in older documentation using \code{permuteGeneral}. This was done for consistency as we felt that the output diverged too much from the general constrained output of \code{permuteGeneral} (See \href{https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#output-order-with-permutegeneral}{Output Order with permuteGeneral}). Research is in the initial stage that is focused on implementing a new family of functions similar to the partition family of functions, but for compositions.
        }
    }
}

\section{Changes in RcppAlgos version 2.4.3 (Release date: 2021-05-30)}{
    \itemize{
        \item Other:
        \itemize{
            \item Fixed old urls
        }
    }
}

\section{Changes in RcppAlgos version 2.4.2 (Release date: Never Released)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added function \code{comboGrid} for efficiently generating the Cartesian product where order does not matter.
        }
        \item Enhancements:
        \itemize{
            \item Refactored code base to reduce binary size
        }
        \item Other:
        \itemize{
            \item Removed LazyData from DESCRIPTION. Also added \code{rmarkdown} to Suggests.
        }
    }
}

\section{Changes in RcppAlgos version 2.4.1 (Release date: 2020-03-24)}{
    \itemize{
        \item Bug Fixes:
        \itemize{
            \item Removed broken URL in NEWS.Rd
        }
    }
}

\section{Changes in RcppAlgos version 2.4.0 (Release date: Never Released)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added \code{comboIter} and \code{permuteIter}. These functions return iterators for iterating over combinations and permutations. They have a similar interface to \code{comboGeneral} and \code{permuteGeneral} and currently only work with standard combinations and permutations. They do not yet work with constraints (This will be the focus of the next release).
            \item Added "High Performance Benchmarks" and "Combinatorial Iterators in RcppAlgos" vignettes
        }
        \item Enhancements:
        \itemize{
            \item Now able to interrupt general constraint problems
        }
        \item Bug Fixes:
        \itemize{
            \item In 2.3.5 and 2.3.6, we mistakingly allowed a \code{constraintFun} to be applied to a logical vector which was crashing R. We have corrected this in 2.4.0+.
            \item Changed the data type for sizing the index matrix in \code{permuteGeneral}. Originally, we were using \code{int} and when the output was large enough, it was causing an integer overflow thus causing the index matrix to be improperly sized. We have since changed the data type to the recommended \code{std::size_t}.
        }
    }
}

\section{Changes in RcppAlgos version 2.3.6 (Release date: 2020-01-29)}{
    \itemize{
        \item Bug Fixes:
        \itemize{
            \item Fixed bug associated with integer vectors, multisets, and constraints. See \href{https://github.com/jwood000/RcppAlgos/issues/12}{Issue #12} for more information.
        }
    }
}

\section{Changes in RcppAlgos version 2.3.5 (Release date: 2020-01-27)}{
    \itemize{
        \item New Features:
        \itemize{
            \item Added \code{comboGroups}, \code{comboGroupsCount}, and \code{comboGroupsSample}. These functions deal with partitioning a vector/set into groups of equal size. See \href{https://stackoverflow.com/q/57732672/4408538}{Create Combinations in R by Groups}. See the related integer sequences A025035-A025042 at \href{https://oeis.org}{OEIS} (E.g. \href{https://oeis.org/A025036}{A025036} for Number of partitions of { 1, 2, ..., 4n } into sets of size 4.)
            \item Added vignettes (First version with vignettes)
            \item Added website via the excellent package pkgdown
            \item Now using C++11 instead of C++14. See \href{https://github.com/jwood000/RcppAlgos/issues/10}{Issue #10} for more information.
        }
        \item Enhancements:
        \itemize{
            \item Extended general partitions algorithm to multisets. E.g. \code{comboGeneral(10, 8, freqs = rep(1:5, 2), constraintFun = "sum", comparisonFun = "==", limitConstraints = 55)}
            \item Improved constraint algorithm for the general case.
            \item Added support for complex and raw types for all combinatorial functions.
            \item Improved permutation algorithm for all cases (I.e. no rep, with rep, and with multisets).
            \item Added loop unrolling to prime sieve algorithm for improved efficiency.
        }
        \item Bug Fixes:
        \itemize{
            \item Corrected checks for total number of partitions and assignment of number of rows when \code{upper} is applied in \code{{combo|permute}General}. See \href{https://github.com/jwood000/RcppAlgos/issues/9#issue-467908452}{Issue #9} for more information.
            \item \code{permuteGeneral} no longer alters source vector. See \href{https://github.com/jwood000/RcppAlgos/issues/11}{Issue #11} for more information.
        }
    }
}

\section{Changes in RcppAlgos version 2.3.4 (Release date: 2019-07-02)}{
    \itemize{
        \item Fixed clang/gcc-ASAN and valgrind issues in 2.3.3. These issues were arising from finding the first vector to meet the criteria in PartitionRep/Distinct. We also found further issues in the standard functions when the length of the partition was 2. The algorithm would eventually try and access an element of a vector at index -1. These fixes were confirmed by successful Rdevel CMD check under 'r-devel-ubsan-clang' via docker using the advice found here: https://knausb.github.io/2017/06/reproducing-a-clang-ubsan-issue/ and http://dirk.eddelbuettel.com/code/sanitizers.html.
    }
}

\section{Changes in RcppAlgos version 2.3.3 (Release date: 2019-06-30)}{
    \itemize{
        \item Fixed clang-UBSAN issue in 2.3.2. It was caused by populating a vector of ints with values larger than 2^31 - 1.
        \item Added optimized algorithm to \code{{combo|permute}General} when constraintFun = "sum", comparisonFun = "==", and the vector passed has a special form. This problem is a special case of the subset sum problem.
        \item Using std::vector and push_back member function instead of pre-allocating matrix when constraint is applied in \code{{combo|permute}General}. This alleviates the need to guess the upper limit and subsequently subset as only elements that meet the constraints are added.
        \item Fixed error in PollardRho.cpp when number passed had factors close to the limit in the predefined lookup table (i.e. constexpr int64_t FirstOmittedPrime = 3989)
    }
}

\section{Changes in RcppAlgos version 2.3.2 (Release date: 2019-04-10)}{
    \itemize{
        \item Fixed clang-UBSAN issue in 2.3.1. It was caused by casting extremely large values to int64_t.
        \item Corrected handling of small values in PrimeSieveBig.
    }
}

\section{Changes in RcppAlgos version 2.3.1 (Release date: 2019-03-20)}{
    \itemize{
        \item Explicitly casted to double for sqrt to silence Solaris
        \item Corrected handling of NaNs
    }
}

\section{Changes in RcppAlgos version 2.3.0 (Release date: 2019-03-17)}{
    \itemize{
        \item All functions now have parallel capabilites via RcppThread.
        \item Utilizes RMatrix.h from RcppParallel for thread safe matrix class.
        \item Major overhaul of primeSieve for large primes.
        \item Added stdThreadMax for obtaining the number of threads available on a machine
    }
}

\section{Changes in RcppAlgos version 2.2.0 (Release date: 2018-09-15)}{
    \itemize{
        \item Disabled \code{Parallel} argument as it was causing unpredictable errors on certain platforms. Research is ongoing to correct this for use in future versions. The development version will retain this feature.
        \item Corrected UBSAN error that caused by filling a vector of integers with numbers larger than 2^31 - 1.
    }
}

\section{Changes in RcppAlgos version 2.1.0 (Release date: 2018-09-11)}{
    \itemize{
        \item Added argument \code{Parallel} to general and sampling functions for increased gains in efficiency.
            \itemize{
                \item \code{comboGeneral(30, 10, Parallel = TRUE)}
                \item \code{permuteGeneral(12, 7, TRUE, constraintFun = "sum", Parallel = TRUE)}
            }
        \item Logical class is now preserved in combinatorial functions
        \item {Added gmp support to combinatorial functions. Now, one can accurately and quickly work with combinations/permutations of large vectors \emph{E.g.}:
            \itemize{
                \item \code{comboSample(runif(10000), 100, n = 10, seed = 42, Parallel = TRUE)}
                \item \code{permuteGeneral(factor(state.name), 20,
                            lower = 1e15, upper = 1e15 + 1000)}
            }
        }
        \item {Added \code{FUN} argument to all combinatorial functions. Allows user to pass custom functions to be applied to combinations/permutations.}
    }
}

\section{Changes in RcppAlgos version 2.0.3 (Release date: 2018-06-30)}{
    \itemize{
        \item Corrected clang UBSAN error identified by two different unit tests. In both situations, the problem was occurring as a result of populating a vector of integers with values from a vector of doubles that contained a nan (Not-a-Number). Most information was obtained from Brian J. Knaus's blog titled : "Reproducing a clang-UBSAN issue" (https://knausb.github.io/2017/06/reproducing-a-clang-ubsan-issue/)
    }
}

\section{Changes in RcppAlgos version 2.0.2 (Release date: 2018-06-15)}{
    \itemize{
        \item Corrected divide by zero in divisorsRcpp unit test.
    }
}

\section{Changes in RcppAlgos version 2.0.1 (Release date: 2018-06-10)}{
    \itemize{
        \item Corrected spelling in DESCRIPTION
    }
}

\section{Changes in RcppAlgos version 2.0.0 (Release date: Never Released)}{
    \itemize{
        \item Changed max value and explicitly casted a few values to (int64_t) in PollardRho.cpp for efficiency while still maintaining accuracy from 2^60 to 2^62 (isPrimeRcpp is roughly 10\% faster now).
        \item Updated core permutation algorithm for greater efficiency and generality
        \item Added capability of generating specific combinations/permutations
        \item Changed arguments to comboGeneral/permuteGeneral. rowCap is now upper.
        \item Added comboSample, permuteSample, comboCount, and permuteCount
        \item Fixed bug with numDivisorSieve and divisorSieve when the lower bound was greater than 1 and less than the sqrt of the upper bound. In the previous version, the numbers in this range would have duplicated values/counts.
        \item Increased efficiency of numDivisorSieve by a factor of 2.
        \item Updated unit tests for greater coverage. See the function package_coverage from the package 'covr'.
    }
}

\section{Changes in RcppAlgos version 1.0.1 (Release date: 2018-03-09)}{
    \itemize{
        \item Corrected precision limits in documentation from 2^64 to 2^63.
        \item Changed const type in PollardRho.cpp from (int64_t) to (double) to correct "UndefinedBehaviorSanitizer"
        \item Changed examples in primeFactorizeSieve to reduce check time
        \item Added RcppAlgos-package man file.
    }
}

\section{Changes in RcppAlgos version 1.0.0 (Release date: 2018-03-07)}{
    \itemize{
        \item Added the following functions: primeFactorize (vectorized pollard rho factorization), divisorsRcpp (vectorized factoring (complete)), isPrimeRcpp (vectorized primality testing using Miller-Rabin algorithm), & primeCount (based on the primecount algorithm by Kim Walisch)
        \item Completely revamped the primeSieve function. It is now a segmented sieve of Eratosthenes with wheel factorization based on primesieve by Kim Walisch.
        \item Renamed divisorsList to divisorsSieve (reason for the major version update to 1.0.0)
        \item Renamed primeFactorizationList to primeFactorizeSieve
        \item Made the sieving functions more flexible. They are now able to generate results over a range and can also produce named objects.
        \item All number theoretic functions have been made more efficient. Some make use of the fast integer division library 'libdivide' by ridiculousfish.
    }
}

\section{Changes in RcppAlgos version 0.2.5 (Release date: 2018-01-04)}{
    \itemize{
        \item Added unit tests.
        \item Removed unnecessary files.
        \item Fixed bug in primeSieve that occurred when a number with a decimal was passed (e.g. 2.01).
        \item Adjusted accepted lower bound for numDivisorSieve, eulerPhiSieve, divisorsList, primeSieve, and primeFactorizationList.
        \item Fixed bug when non-unique elements are present with factors.
    }
}

\section{Changes in RcppAlgos version 0.2.4 (Release date: 2017-12-18)}{
    \itemize{
        \item Fixed bug that occurs when non-unique elements are present for combinations with replacement.
    }
}

\section{Changes in RcppAlgos version 0.2.3 (Release date: 2017-12-18)}{
    \itemize{
        \item Fixed segmentation fault error highlighted by valgrind check in version 0.2.2.
        \item Updated DESCRIPTION file.
    }
}

\section{Changes in RcppAlgos version 0.2.2 (Release date: 2017-12-15)}{
    \itemize{
        \item Fixed bug in constraint functions that occurred when m = 1 and the constraint limit was equal to the last element in v. It was returning a 2x1 matrix with the same value twice.  It is now correctly returning a 1x1 matrix with the correct value 1 time.
        \item Reorganized source code such that all utility functions for the combinatoric functions are now in their own file. Additionally added header for this file.
        \item All combinatoric functions can now utilize the rowCap argument. Before, rowCap only applied to the combinatorial functions with constraints.
        \item comboGeneral can now find all combinations of multisets.
        \item Both comboGeneral and permuteGeneral can utilize the argument m when dealing with multisets.  Before, permuteGeneral would simply return all permutations of a multiset. Now you can specify the lengths of the output.
    }
}

\section{Changes in RcppAlgos version 0.2.1 (Release date: 2017-11-29)}{
    \itemize{
        \item Fixed bug that would occur in two edge cases involving the constraint functions.
        \item Slightly modified formatting for primeSieve.Rd
    }
}

\section{Changes in RcppAlgos version 0.2.0 (Release date: 2017-11-28)}{
    \itemize{
        \item Updated combination algorithms. They are now more than twice as fast.
        \item Updated constraint functions so that memory access is always within container bounds
        \item Consolidated redundant code for outputting different Rcpp types (e.g. IntegerMatrix, CharacterMatrix, etc.) via a templated approach.
        \item Added the function permuteGeneral that is analogous to comboGeneral only instead of combinations, it gives all permutations. It has an additional argument (i.e. 'freqs') that is used to generate permutations of multisets.
        \item All combinatoric functions now support factor types.
    }
}

\section{Changes in RcppAlgos version 0.1.2 (Release date: 2017-11-03)}{
    \itemize{
        \item Corrected minor typo in README file.
        \item Fixed minor error regarding explicitly comparing variables to large numbers that are typed out. Simply adding a decimal along with a zero remedies the situation.
    }
}

\section{Changes in RcppAlgos version 0.1.1 (Release date: 2017-11-03)}{
    \itemize{
        \item Improved ComboConstraint function by removing unnecessary subsetting.
        \item Improved PrimeSieve internal C++ algorithm.
        \item Corrected the errors with respect to the math functions in C++. Explicitly overloaded the parameters of these functions by casting them to the double type.
    }
}

\section{RcppAlgos version 0.1.0 (Release date: 2017-10-26)}{
    \itemize{
        \item Initial Release
    }
}
