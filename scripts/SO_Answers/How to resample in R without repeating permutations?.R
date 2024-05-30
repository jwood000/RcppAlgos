reprex::reprex({
    #' The package `RcppAlgos` (I am the author) has some functions specifically for this type of problem. Observe:

    library(RcppAlgos)
    permuteSample(v = 0:10, m = 5, seed = 123, n = 3)

    ## Or use global RNG to get the same results
    set.seed(123)
    permuteSample(v = 0:10, m = 5, n = 3)

    #' They are guaranteed to sample **without** replacement:

    permuteCount(3)
    permuteSample(3, n = 7)

    #' If we request the maximal number of samples, it will be equivalent to generating all results and shuffling:

    ps <- permuteSample(3, n = 6, seed = 100)

    ## Shuffled
    ps

    ## Same as getting them all
    identical(
        permuteGeneral(3),
        ps[do.call(order, as.data.frame(ps)), ]
    )

    #' ## Flexibility
    #'
    #' These functions are quite general as well. For example, if we want to sample a vector where repetition _within the vector_ is allowed set `repetition = TRUE` (not to be confused with sampling with replacement):

    permuteSample(v = 0:10, m = 5, repetition = TRUE, seed = 42, n = 3)

    ## Go beyond the vector length
    permuteSample(v = 3, m = 10, repetition = TRUE, seed = 42, n = 3)

    #' What about specific multiplicity for each element... easy! Use the `freqs` parameter:

    ## With v = 3 and freqs = 2:4, this means:
    ##
    ##    1 can occur a maximum of 2 times
    ##    2 can occur a maximum of 3 times
    ##    3 can occur a maximum of 4 times
    ##
    permuteSample(v = 3, m = 6, freqs = 2:4, seed = 1234, n = 3)

    #' It works on all atomic types as well:

    ## When m is NULL, the length defaults to the length of v
    set.seed(28)
    permuteSample(
        v = rnorm(3) + rnorm(3) * 1i,
        repetition = TRUE,
        n = 5
    )

    #' Sampling With Replacement
    #'
    #' The default behavior of the sampling functions in `RcppAlgos` is to sample without replacement. If you need sampling with replacement, we can make use of the `sampleVec` parameter. This parameter takes a vector and interprets them as indices representing the lexicographical permutations to return. For example:

    permuteSample(3, sampleVec = c(1, 3, 6))

    ## Same as
    permuteGeneral(3)[c(1, 3, 6), ]

    #' This means we can generate our own random sample of indices and pass it to `permuteSample`:

    set.seed(123)
    idx <- sample(permuteCount(3), 5, replace = TRUE)
    idx

    permuteSample(3, sampleVec = idx)

    ## See rownames corresponding to the lexicographical result
    permuteSample(3, sampleVec = idx, namedSample = TRUE)

    #' ## Efficiency
    #'
    #' The functions are written in `C++` with efficiency in mind. Take for example the largest case @whuber tests:

    system.time(size_10 <- permuteSample(10, n = 1e6, seed = 42))

    ## Guaranteed no repeated permutations!
    nrow(unique(size_10))

    system.time(size_15 <- permuteSample(15, n = 1e6, seed = 42))

    ## Guaranteed no repeated permutations!
    nrow(unique(size_15))

    #' Why not try even larger cases?!?!?! Using `nThreads`, we can achieve even greater efficiency:

    ## Single threaded
    system.time(permuteSample(15, n = 1e7, seed = 321))

    ## Using 4 threads
    system.time(permuteSample(15, n = 1e7, seed = 321, nThreads = 4))

    #' ## Sampling Other Combinatorial Objects
    #'
    #' There are also sampling functions for other combinatorial objects:
    #'
    #'   1. permutations: `permuteSample`
    #'   2. combinations: `comboSample`
    #'   3. integer partitions: `partitionsSample`
    #'   4. integer compositions: `compositionsSample`
    #'   5. partitions of a vector into groups `comboGroupsSample`

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
