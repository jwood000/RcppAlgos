reprex::reprex({
    #'
    #' This document covers the topic of [integer compositions](<https://en.wikipedia.org/wiki/Composition_(combinatorics)>) in `RcppAlgos`. Integer compositions are closely related to integer partitions; however, order matters. While standard compositions are well studied and widely implemented, compositions with distinct parts present additional structural constraints that require specialized treatment. To our knowledge, `RcppAlgos` is the first combinatorics library to provide a dedicated next-lexicographical algorithm for efficiently generating distinct integer compositions, enabling fast enumeration, ranking, and sampling even for large problems.
    #'
    #' ***
    #'
    #' # Integer Compositions
    #'
    #' Compositions are generated using `compositionsGeneral()`, which follows the same interface as the other combinatorial functions in `RcppAlgos`.
    #'
    #' Prior to version `2.10.0`, the composition engine in `RcppAlgos` supported only a restricted subset of compositions with repetition.
    #'
    #' Version `2.10.0` introduces complete support for compositions with repetition and adds full support for compositions with distinct parts.
    #'
    #' Planned future enhancements include comprehensive support for compositions of multisets.
    #'
    #' The output with `compositionsGeneral` is in lexicographical order across
    #' all matching compositions. Its `lower` and `upper` index those matching
    #' rows. When matches exist, `lower = 1` gives the same prefix as
    #' omitting `lower`. With no matches, omitting both bounds returns an
    #' empty matrix; explicit positive bounds exceed the zero result count
    #' and raise an error.
    #' When `weak = TRUE` and zero is included, zeros are actual parts of the
    #' sequence: `c(0, 0, 5)`, `c(0, 5, 0)`, and `c(5, 0, 0)` are weak
    #' compositions of 5. Increasing the width permits additional zeros.
    #'
    #' For more general cases, we can use `permuteGeneral`. For supported
    #' sum/equality cases without `lower`, it uses partition-then-permute order,
    #' even when an equivalent composition algorithm exists. Supplying `lower`
    #' instead filters ordinary permutation ranks. The section below explains
    #' these contracts. In permutations, zero is an actual element that can
    #' occupy any position, subject to its multiplicity.
    #'
    #' With that in mind, generating compositions with `RcppAlgos` is easy, flexible, and quite efficient.
    #'
    #' We continue to use the `ht` function defined in the [Combination and Permutation Basics](<https://jwood000.github.io/RcppAlgos/articles/GeneralCombinatorics.html>) vignette):

    ht <- function(d, m = 5, n = m) {
        ## print the head and tail together
        cat("head -->\n")
        print(head(d, m))
        cat("--------\n")
        cat("tail -->\n")
        print(tail(d, n))
    }

    #'
    #' ## Standard Compositions
    #'
    #' As with partitions, the standard definition of a composition allows parts to repeat and represents the most commonly encountered form.
    #'
    #' In this section we will explore compositions of this form.
    #'
    #' ### Case 1: All Compositions of _N_
    #'

    library(RcppAlgos)
    options(width = 90)

    packageVersion("RcppAlgos")

    cat(paste(capture.output(sessionInfo())[1:3], collapse = "\n"))

    compositionsGeneral(0:3, repetition = TRUE)

    ## Weak compositions
    compositionsGeneral(0:3, repetition = TRUE, weak = TRUE)

    ## Get weak compositions with width > than target
    ht(compositionsGeneral(0:3, 10, repetition = TRUE, weak = TRUE))

    #'
    #' ### Case 2: Compositions of _N_ of Length _m_
    #'

    compositionsGeneral(6, 3, repetition = TRUE)

    ## Including zero
    ht(compositionsGeneral(0:6, 3, repetition = TRUE))

    ## Weak compositions of length 3
    ht(compositionsGeneral(0:6, 3, repetition = TRUE, weak = TRUE))

    #'
    #' ## Distinct Compositions
    #'
    #' Distinct compositions impose the additional constraint that all parts must be unique. While conceptually simple, this restriction substantially complicates efficient generation, counting, ranking, and sampling.
    #'
    #' Version `2.10.0` introduces a next-lexicographical generation algorithm for distinct integer compositions, enabling high-performance enumeration and large-scale computations that were previously computationally prohibitive. To our knowledge, this represents the first implementation of a fully integrated next-lex engine for distinct compositions in an open-source combinatorics library.
    #'
    #' ### Case 3: Compositions of _N_ into Distinct Parts
    #'

    compositionsGeneral(6)

    ht(compositionsGeneral(40, 7))

    ## Note the subtlety here. When zero is NOT included, the last result is
    ## simply the reverse of the first result. Below, this is not the case.
    ht(compositionsGeneral(0:40, 7))

    #'
    #' #### Verifying Distinct Compositions via Partitions (Brute Force Construction)
    #'
    #' This section verifies correctness and illustrates the combinatorial overhead involved in generating distinct compositions in lexicographical order using a brute-force construction.
    #'
    #' We will use the example above of generating distinct compositions of 40 into 7 parts (including zero, non-weak case).
    #'
    #' Since we are not working with weak compositions, zero cannot appear as an ordinary part. To account for zero-padding, we proceed as follows:
    #'
    #' * Generate all distinct partitions of 40 of length 6.
    #' * Permute each partition to obtain the corresponding compositions.
    #' * `cbind` a zero column.
    #' * Repeat for partitions of length 7 (no zero padding required).
    #' * Combine and lexicographically sort the results.
    #'

    ## Helper: permute each partition to obtain all corresponding compositions
    get_perms_of_parts <- function(parts) {
        do.call(
            rbind,
            lapply(seq_len(nrow(parts)), \(i) permuteGeneral(table(parts[i, ])))
        )
    }

    ## Length-6 partitions become length-7 compositions by prepending one zero
    perms_6 <- get_perms_of_parts(partitionsGeneral(40, 6))
    res_6   <- cbind(0L, perms_6)

    ## Length-7 partitions already match the desired output width
    res_7 <- get_perms_of_parts(partitionsGeneral(40, 7))

    ## Combine and sort in lexicographical order
    res       <- rbind(res_6, res_7)
    brute_lex <- res[do.call(order, as.data.frame(res)), ]

    ## Verify against the direct generator
    identical(compositionsGeneral(0:40, 7), brute_lex)

    #'
    #' ### Case 4: Integer Compositions of _N_ into Parts of Varying Multiplicity
    #'
    #' Currently, we must use `permuteGeneral`.
    #'

    ## Generates error
    compositionsGeneral(5, 3, freqs = rep(2, 5))

    ## compositions of 5 into 3 parts where each part can
    ## be used a maximum of 2 times.
    permuteGeneral(5, 3, freqs = rep(2, 5),
                   constraintFun = "sum",
                   comparisonFun = "==",
                   limitConstraints = 5)

    #'
    #' ## Generating Compositions with `permuteGeneral()`
    #'
    #' For supported `constraintFun = "sum"` and `comparisonFun = "=="` calls
    #' without `lower`, `permuteGeneral` uses **partition-then-permute order**.
    #' It visits matching partitions in lexicographical order and emits all
    #' distinct permutations of each partition in lexicographical order before
    #' moving to the next. The combined output need not be globally sorted.
    #' These cases include integer-valued, evenly spaced source vectors and an
    #' integer target, with supported repetition, distinct, or multiset rules.
    #'
    #' `upper` alone caps matching rows in that order and may stop partway
    #' through a partition's permutations. If `lower` is supplied, both bounds
    #' instead refer to ordinary permutation ranks before applying the
    #' constraint. The surviving rows are in lexicographical order. Omitting
    #' `upper` tests all candidates from `lower` onward.
    #'
    #' Consequently, `lower = 1` is not redundant. Requesting four matches
    #' with `upper` alone differs from testing permutation ranks 1 through 4.
    #' Neither interpretation makes these bounds composition ranks; use
    #' `compositionsGeneral` to index lexicographically ordered compositions.
    #'

    ## No lower: permute each matching partition before moving to the next.
    permuteGeneral(6, 3, repetition = TRUE, constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 6, upper = 4)

    ## Explicit lower: test ordinary permutation ranks 1 through 4.
    ## Only (1, 1, 4) matches; upper is not a request for four matches.
    permuteGeneral(6, 3, repetition = TRUE, constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 6,
                   lower = 1, upper = 4)

    ## Compositions: return the first four matching rows in lexicographic order.
    compositionsGeneral(6, 3, repetition = TRUE, target = 6, upper = 4)

    #'
    #' Zero is an actual element in `permuteGeneral`, so it may occur anywhere
    #' in a result, subject to multiplicities. This differs from non-weak
    #' composition padding.
    #'
    #' The following calls omit `lower` and therefore retain
    #' partition-then-permute order:
    #'

    ## With repetition
    permuteGeneral(0:3, repetition = TRUE,
                   constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 3)

    ## Similar behavior as with compositionsGeneral when width > target and
    ## weak = TRUE. That is, we generate more results b/c zero is included
    ## in the output sequences
    tail(permuteGeneral(0:3, 10, repetition = TRUE,
                        constraintFun = "sum",
                        comparisonFun = "==", limitConstraints = 3))

    ## Distinct Parts
    permuteGeneral(0:3,
                   constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 3)

    ## Distinct Parts and Specific Size
    permuteGeneral(0:3, 3, repetition = FALSE,
                   constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 3)

    #'
    #' ## The Role of `target`
    #'
    #' Just as with integer partitions, we can make use of the `target` parameter to specify the integer being partitioned.
    #'
    #' Below, we replicate the examples from the [Integer Partitions](<https://jwood000.github.io/RcppAlgos/articles/IntegerPartitions.html#the-role-of-target>) vignette, this time generating compositions. Because compositions yield substantially more results, we use `ht()` to limit the displayed output.
    #'

    ## Here we partition 30 using only distinct parts up to 10
    ht(compositionsGeneral(10, 5, target = 30))

    ## Here we partition 15 using only parts up to 4 with zero included
    ht(compositionsGeneral(0:4, 5, repetition = TRUE, target = 15))

    ## Here we partition 22 using only parts up to 8 with zero(s) included
    ht(compositionsGeneral(0:8, 6, freqs = c(4, rep(1, 8)), target = 22))

    ## Same as above, just making use of the table method
    ht(compositionsGeneral(table(c(rep(0L, 4), 1:8)), 6, target = 22))

    #'
    #' ## Efficiency Generating Partitions and Compositions
    #'
    #' Supported `compositionsGeneral` calls can use parallel generation.
    #' Partition-then-permute generation through `permuteGeneral` without
    #' `lower` is serial. Supplying `lower` allows parallel candidate filtering
    #' when the requested range is large enough; this scans ordinary
    #' permutations and can be expensive when matches are sparse.
    #'

    ## compositions of 25
    system.time(compositionsGeneral(0:25, repetition = TRUE))

    compositionsCount(0:25, repetition = TRUE)

    ## Use multiple threads for greater efficiency. Generate
    ## over 16 million compositions instantly
    system.time(compositionsGeneral(0:25, repetition = TRUE, nThreads = 4))


    ## weak compositions of 12 using nThreads = 4
    system.time(weakComp12 <- compositionsGeneral(0:12, repetition = TRUE,
                                                  weak = TRUE, nThreads = 4))

    ## And using permuteGeneral
    system.time(weakPerm12 <- permuteGeneral(0:12, 12, repetition = TRUE,
                                             constraintFun = "sum",
                                             comparisonFun = "==",
                                             limitConstraints = 12))

    dim(weakPerm12)

    identical(weakPerm12[do.call(order, as.data.frame(weakPerm12)), ],
              weakComp12)

    ## Distinct Parts
    system.time(compositionsGeneral(55, 8))

    compositionsCount(55, 8)

    ## Similar to above... using 4 threads gets us the result instantly
    system.time(compositionsGeneral(55, 8, nThreads = 4))

    ## General compositions with varying multiplicities
    system.time(comp25_gen <- permuteGeneral(25, 10, freqs = rep(4:8, 5),
                                             constraintFun = "sum",
                                             comparisonFun = "==",
                                             limitConstraints = 25))

    dim(comp25_gen)

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
