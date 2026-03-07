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
    #' The output with `compositionGeneral` will be in lexicographical order. When we set `weak = TRUE` and zero is included, we will obtain **_weak compositions_**, which allow for zeros to be a part of the sequence (E.g. `c(0, 0, 5), c(0, 5, 0), c(5, 0, 0)` are weak compositions of 5). As the Wikipedia article points out, we can increase the number of zeros indefinitely when `weak = TRUE`.
    #'
    #' For more general cases, we can make use of `permuteGeneral`, keeping in mind that the output will not be in lexicographical order. Another consideration with `permuteGeneral` is that when we include zero, we will always obtain weak compositions.
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
    #' As noted in the introduction, compositions can also be generated using `permuteGeneral()`. As with the other general combinatorial generators in `RcppAlgos`, the procedure first finds feasible combinations satisfying the constraint and then produces permutations of those combinations. In this setting, the combinations correspond to integer partitions of the target (`limitConstraints` in this case).
    #'
    #' This also implies that we always produce **weak compositions** when zero is included. Furthermore, since the results are not generated in lexicographical order, parameters such as `lower` and `upper` cannot be applied, as there is no well-defined notion of the _n<sup>th</sup>_ result.
    #'
    #' For example, the following generates all weak compositions of 3 using elements from `0:3` and NOT in lex-order:
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
    #' With `compositionGeneral` we are able to take advantage of parallel computation. With `permuteGeneral`, the parallel options have no effect when generating compositions.
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
