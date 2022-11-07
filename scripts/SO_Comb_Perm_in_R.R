reprex::reprex({
    #' # A Walk Through a Slice of Combinatorics in R*
    #'
    #' Below, we examine packages equipped with the capabilities of generating combinations & permutations. If I have left out any package, please forgive me and please leave a comment or better yet, edit this post.
    #'
    #' Outline of analysis:
    #'
    #'   1. Introduction
    #'   2. Setup
    #'   3. Combinations
    #'   4. Permutations
    #'   5. Multisets
    #'   6. Summary
    #'   7. Memory
    #'
    #' Before we begin, we note that combinations/permutations **with** replacement of distinct vs. non-distint items chosen _m_ at a time are equivalent. This is so, because when we have replacement, it is not specific. Thus, no matter how many times a particular element originally occurs, the output will have an instance(s) of that element repeated 1 to _m_ times.
    #'
    #' ## 1. Introduction
    #'
    #' ### Packages:
    #'
    #'   1. `gtools`
    #'   2. `combinat`
    #'   3. `multicool`
    #'   4. `partitions`
    #'   5. `RcppAlgos`
    #'   6. `arrangements`
    #'   7. `utils`
    #'
    #' I did not include `permute` or `permutations` as they are not really meant to attack these types of problems. I also did not include the updated `gRbase` as certain cases crashed my computer.
    #'
    #'
    #' #### |--------------------------------------- **OVERVIEW** ----------------------------------------|
    #'
    #'
    #'     |_________________| gtools | combinat | multicool | partitions |
    #'     |       comb rep  |  Yes   |          |           |            |
    #'     |    comb NO rep  |  Yes   |   Yes    |           |            |
    #'     |       perm rep  |  Yes   |          |           |            |
    #'     |    perm NO rep  |  Yes   |   Yes    |    Yes    |    Yes     |
    #'     |  perm multiset  |        |          |    Yes    |    Yes     |
    #'     |  comb multiset  |        |          |           |            |
    #'     | accepts factors |        |   Yes    |           |            |
    #'     |    m at a time  |  Yes   |  Yes/No  |           |            |
    #'     | general vector  |  Yes   |   Yes    |    Yes    |            |
    #'     |     iterable    |        |          |    Yes    |            |
    #'     | parallelizable  |        |          |           |            |
    #'     | multi-threaded  |        |          |           |            |
    #'     |   big integer   |        |          |           |            |
    #'
    #'     |_________________| arrangements | RcppAlgos | utils |
    #'     |       comb rep  |     Yes      |    Yes    |       |
    #'     |    comb NO rep  |     Yes      |    Yes    |  Yes  |
    #'     |       perm rep  |     Yes      |    Yes    |       |
    #'     |    perm NO rep  |     Yes      |    Yes    |       |
    #'     |  perm multiset  |     Yes      |    Yes    |       |
    #'     |  comb multiset  |     Yes      |    Yes    |       |
    #'     | accepts factors |     Yes      |    Yes    |  Yes  |
    #'     |    m at a time  |     Yes      |    Yes    |  Yes  |
    #'     | general vector  |     Yes      |    Yes    |  Yes  |
    #'     |     iterable    |     Yes      |    Yes    |       |
    #'     | parallelizable  |     Yes      |    Yes    |       |
    #'     |   big integer   |     Yes      |    Yes    |       |
    #'     | multi-threaded  |              |    Yes    |       |
    #'
    #' The tasks, `m at a time` and `general vector`, refer to the capability of generating results "_m_ at a time" and rearranging a "general vector" as opposed to `1:n`. In practice, we are generally concerned with finding rearrangements of a general vector, therefore all examinations below will reflect this when possible.
    #'
    #' ## 2. Setup
    #'
    #' All benchmarks were ran on 3 different set-ups.
    #'
    #'   1. 2022 Macbook Air Apple M2 24 GB
    #'   2. 2020 Macbook Pro i7 16 GB
    #'   3. 2022 Windows Surface i5 16 GB

    library(microbenchmark)
    ## print up to 4 digits to keep microbenchmark output tidy
    options(digits = 4)
    options(width = 90)

    numThreads <- min(as.integer(RcppAlgos::stdThreadMax() / 2), 6)
    numThreads

    pkgs <- c("gtools", "combinat", "multicool", "partitions",
              "RcppAlgos", "arrangements", "utils", "microbenchmark")
    sapply(pkgs, packageVersion, simplify = FALSE)

    #' The listed results were obtained from setup #1 (i.e. Macbook Air M2). The results on the Macbook Pro were similar, however with the Windows setup, multi-threading was less effective. In some cases on the Windows setup, the serial execution was faster. We will call all functions with the pattern `package::function` so no `library` calls are needed.
    #'
    #' ## 3. Combinations
    #'
    #' First, we examine combinations without replacement chosen _m_ at a time.
    #'
    #'   1. `RcppAlgos`
    #'   2. `combinat`
    #'   3. `gtools`
    #'   4. `arrangements`
    #'   5. `utils`
    #'
    #' How to:

    set.seed(13)
    tVec1 <- sort(sample(100, 20))
    m <- 10
    t1 <- RcppAlgos::comboGeneral(tVec1, m)  ## returns matrix with m columns
    t3 <- combinat::combn(tVec1, m)  ## returns matrix with m rows
    t4 <- gtools::combinations(20, m, tVec1)  ## returns matrix with m columns
    identical(t(t3), t4) ## must transpose to compare
    t5 <- arrangements::combinations(tVec1, m)
    identical(t1, t5)
    t6 <- utils::combn(tVec1, m)  ## returns matrix with m rows
    identical(t(t6), t4) ## must transpose to compare

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::comboGeneral(tVec1, m,
                                                 nThreads = numThreads),
        cbRcppAlgosSer = RcppAlgos::comboGeneral(tVec1, m),
        cbGtools = gtools::combinations(20, m, tVec1),
        cbCombinat = combinat::combn(tVec1, m),
        cbUtils = utils::combn(tVec1, m),
        cbArrangements = arrangements::combinations(tVec1, m),
        unit = "relative"
    )

    #' Now, we examine combinations with replacement chosen _m_ at a time.
    #'
    #'   1. `RcppAlgos`
    #'   2. `gtools`
    #'   3. `arrangements`
    #'
    #' How to:

    set.seed(97)
    tVec2 <- sort(rnorm(12))
    m  <- 9
    t1 <- RcppAlgos::comboGeneral(tVec2, m, repetition = TRUE)
    t3 <- gtools::combinations(12, m, tVec2, repeats.allowed = TRUE)
    identical(t1, t3)
    t4 <- arrangements::combinations(tVec2, m, replace = TRUE)
    identical(t1, t4)

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::comboGeneral(tVec2, m, TRUE,
                                                 nThreads = numThreads),
        cbRcppAlgosSer = RcppAlgos::comboGeneral(tVec2, m, TRUE),
        cbGtools = gtools::combinations(12, m, tVec2, repeats.allowed = TRUE),
        cbArrangements = arrangements::combinations(tVec2, m, replace = TRUE),
        unit = "relative"
    )

    #' ## 4. Permutations
    #' First, we examine permutations without replacement chosen _m_ at a time.
    #'
    #'   1. `RcppAlgos`
    #'   2. `gtools`
    #'   3. `arrangements`
    #'
    #' How to:

    tVec3 <- as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29))
    t1 <- RcppAlgos::permuteGeneral(tVec3, 6)
    t3 <- gtools::permutations(10, 6, tVec3)
    identical(t1, t3)
    t4 <- arrangements::permutations(tVec3, 6)
    identical(t1, t4)

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(tVec3, 6,
                                                   nThreads = numThreads),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(tVec3, 6),
        cbGtools = gtools::permutations(10, 6, tVec3),
        cbArrangements = arrangements::permutations(tVec3, 6),
        unit = "relative"
    )

    #' Next, we examine permutations without replacement with a general vector (returning all permutations).
    #'
    #'   1. `RcppAlgos`
    #'   2. `gtools`
    #'   3. `combinat`
    #'   4. `multicool`
    #'   5. `arrangements`
    #'
    #' How to:

    tVec3Prime <- tVec3[1:9]
    ## For RcppAlgos, arrangements, & gtools (see above)

    t4 <- combinat::permn(tVec3Prime) ## returns a list of vectors
    ## convert to a matrix
    t4 <- do.call(rbind, t4)
    t5 <- multicool::allPerm(multicool::initMC(tVec3Prime)) ## returns a matrix with n columns
    all.equal(t4[do.call(order,as.data.frame(t4)),],
              t5[do.call(order,as.data.frame(t5)),])

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(tVec3Prime,
                                                   nThreads = numThreads),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(tVec3Prime),
        cbGtools = gtools::permutations(9, 9, tVec3Prime),
        cbCombinat = combinat::permn(tVec3Prime),
        cbMulticool = multicool::allPerm(multicool::initMC(tVec3Prime)),
        cbArrangements = arrangements::permutations(tVec3Prime),
        times = 25,
        unit = "relative"
    )

    #' Now, we examine permutations without replacement for `1:n` (returning all permutations).
    #'
    #'   1. `RcppAlgos`
    #'   2. `gtools`
    #'   3. `combinat`
    #'   4. `multicool`
    #'   5. `partitions`
    #'   6. `arrangements`
    #'
    #' How to:

    t1 <- partitions::perms(9)  ## returns an object of type 'partition' with n rows
    identical(t(as.matrix(t1)), RcppAlgos::permuteGeneral(9))

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(9, nThreads = numThreads),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(9),
        cbGtools = gtools::permutations(9, 9),
        cbCombinat = combinat::permn(9),
        cbMulticool = multicool::allPerm(multicool::initMC(1:9)),
        cbPartitions = partitions::perms(9),
        cbArrangements = arrangements::permutations(9),
        times = 25,
        unit = "relative"
    )

    #' Lastly, we examine permutations with replacement.
    #'
    #'   1. `RcppAlgos`
    #'   2. `gtools`
    #'   3. `arrangements`
    #'
    #' How to:

    t1 <- RcppAlgos::permuteGeneral(tVec3, 5, repetition = TRUE)
    t3 <- gtools::permutations(10, 5, tVec3, repeats.allowed = TRUE)
    t4 <- arrangements::permutations(x = tVec3, k = 5, replace = TRUE)

    identical(t1, t3)
    identical(t1, t4)

    #' This next benchmark is a little surprising given the results up until now.

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(
            tVec3, 5, TRUE, nThreads = numThreads
        ),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(tVec3, 5, TRUE),
        cbGtools = gtools::permutations(10, 5, tVec3, repeats.allowed = TRUE),
        cbArrangements = arrangements::permutations(tVec3, 5, replace = TRUE),
        unit = "relative"
    )

    #' That is not a typo... `gtools::permutations` is almost as fast as the other compiled functions. I encourage the reader to go check out the source code of `gtools::permutations` as it is one of the most elegant displays of programming out there (`R` or otherwise).
    #'
    #' ## 5. Multisets
    #'
    #' First, we examine combinations of multisets.
    #'
    #'   1. `RcppAlgos`
    #'   2. `arrangements`
    #'
    #' To find combinations/permutations of multisets, with `RcppAlgos` use the `freqs` arguments to specify how many times each element of the source vector, `v`, is repeated.

    set.seed(496)
    myFreqs <- sample(1:5, 12, replace = TRUE)
    ## This is how many times each element will be repeated
    tVec4 <- as.integer(c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233))
    t1 <- RcppAlgos::comboGeneral(tVec4, 12, freqs = myFreqs)
    t3 <- arrangements::combinations(tVec4, 12, freq = myFreqs)
    identical(t1, t3)

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::comboGeneral(
            tVec4, 12, freqs = myFreqs, nThreads = numThreads
        ),
        cbRcppAlgosSer = RcppAlgos::comboGeneral(tVec4, 12, freqs = myFreqs),
        cbArrangements = arrangements::combinations(tVec4, 12, freq = myFreqs),
        unit = "relative"
    )

    #' For permutations of multisets chosen _m_ at a time, we have:
    #'
    #'   1. `RcppAlgos`
    #'   2. `arrangements`
    #'
    #' How to:

    set.seed(123)
    tVec5 <- sort(runif(5))
    t1 <- RcppAlgos::permuteGeneral(tVec5, 8, freqs = rep(4, 5))
    t3 <- arrangements::permutations(tVec5, 8, freq = rep(4, 5))
    identical(t1, t3)

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(
            tVec5, 8, freqs = rep(4, 5), nThreads = numThreads
        ),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(tVec5, 8, freqs = rep(4, 5)),
        cbArrangements = arrangements::permutations(tVec5, 8, freq = rep(4, 5)),
        unit = "relative"
    )

    #' For permutations of multisets returning all permutations, we have:
    #'
    #'   1. `RcppAlgos`
    #'   2. `multicool`
    #'   3. `partitions`
    #'   4. `arrangements`
    #'
    #' How to:

    tVec6 <- (1:5)^3
    ## For multicool, you must have the elements explicitly repeated
    tVec6Prime <- rep(tVec6, times = rep(2, 5))

    ## for comparison
    t1 <- RcppAlgos::permuteGeneral(tVec6, freqs = rep(2, 5))
    t2 <- partitions::multiset(tVec6Prime)
    t3 <- multicool::allPerm(multicool::initMC(tVec6Prime))
    t4 <- arrangements::permutations(tVec6, freq = rep(2, 5))

    ## the package partitions, returns class of integer
    ## whereas RcppAlgos preserves class of tVec6 (i.e. numeric)
    all.equal(t1, t(as.matrix(t2)))
    identical(t1[do.call(order,as.data.frame(t1)),],
              t3[do.call(order,as.data.frame(t3)),])
    identical(t1, t4)

    #' Benchmark:

    microbenchmark(
        cbRcppAlgosPar = RcppAlgos::permuteGeneral(
            tVec6, freqs = rep(2, 5), nThreads = numThreads
        ),
        cbRcppAlgosSer = RcppAlgos::permuteGeneral(tVec6, freqs = rep(2, 5)),
        cbMulticool = multicool::allPerm(multicool::initMC(tVec6Prime)),
        cbPartitions = partitions::multiset(tVec6Prime),
        cbArrangements = arrangements::permutations(tVec6, freq = rep(2, 5)),
        unit = "relative"
    )

    #' ## 6. Summary
    #'
    #' Both `gtools` and `combinat` are well established packages for rearranging elements of a vector. With `gtools` there are a few more options (see the overview above) and with `combinat`, you can rearrange `factors`. With `multicool`, one is able to rearrange multisets. Although `partitions` is limited for the purposes of this question, it is a powerhouse packed with highly efficient functions for dealing with integer partitions.
    #'
    #' `arrangements`
    #' --------
    #'
    #'   1. The output is in lexicographical order.
    #'   2. Allows the user to specify the format via the `layout` argument ("row : row-major", "colmnn : column-major", and "list : list").
    #'   3. Offers convenient methods such as `collect` & `getnext` when working with iterators.
    #'   4. Allows for the generation of more than `2^31 - 1` combinations/permutations via `getnext`. N.B. `RcppAlgos` (via `nextItem`) and `multicool` (via `nextPerm`) are also capable of doing this.
    #'   5. GMP support allows for exploration of combinations/permutations of vectors with many results.
    #'
    #' Observe:

    icomb <- arrangements::icombinations(1000, 7)
    icomb$getnext()

    icomb$getnext(d = 5)

    #' This feature is really nice when you only want a few combinations/permutations. With traditional methods, you would have to generate all combinations/permutations and then subset. This would render the previous example impossible as there are more than `10^17` results (i.e. `ncombinations(1000, 7, bigz = TRUE)` = 194280608456793000).
    #'
    #' This feature along with the improvements to the generators in `arrangements`, allow it to be very efficient with respect to memory.
    #'
    #' `RcppAlgos`
    #' ----------
    #'
    #'   1. The output is in lexicographical order.
    #'   2. There are convenient constraint features that we will not discuss here as they are off-topic for this question. I will only note that the types of problems that can be solved by utilizing these features were the motivation for creating this package (partitions, subset-sum, etc.).
    #'   3. GMP support allows for exploration of combinations/permutations of vectors with many results.
    #'   4. Produce results in parallel using the `Parallel` or `nThreads` arguments.
    #'   5. Similar to `combn`, there is a `FUN` argument for applying a function to each result (See also `FUN.VALUE`).
    #'   6. Provides flexible and merory efficient iterators that allow for bidirectional iteration as well as random access.
    #'      - `nextItem`|`nextNIter`|`nextRemaining`: Retrieve the _next_ lexicographical result(s)
    #'      - `prevItem`|`prevNIter`|`prevRemaining`: Retrieve the _previous_ lexicographical result(s)
    #'      - `front`|`back`|`[[`: Random access methods
    #'      - Allows for easy generation of more than `2^31 - 1` results from any starting place.
    #'
    #' Observe:

    iter <- RcppAlgos::comboIter(1000, 7)

    # first combinations
    iter@nextIter()

    # next 5 combinations
    iter@nextNIter(5)

    # from the current state, the previous combination
    iter@prevIter()

    # the last combination
    iter@back()

    # the 5th combination
    iter[[5]]

    # you can even pass a vector of indices
    iter[[c(1, 3, 5)]]

    # start iterating from any index
    iter[[gmp::pow.bigz(2, 31)]]

    # get useful info about the current state
    iter@summary()

    ## get next ieteration
    iter@nextIter()

    #' In case you were wondering how each package scales, I will leave you with this final example that measures how fast `RcppAlgos` and the `arrangements` packages can generate over 100 million results. Note, `gtools::combinations` is left out here as it will throw the error: `evaluation nested too deeply...`. We also leave out `combn` as it takes quite some time to execute. Curiously, the differences in memory usage between `utils::combn` and `combinat::combn` is quite bizarre given that they are only marginally different (see `?utils::combn` under the "Authors" section).
    #'
    #' Observe:

    set.seed(2187)
    tVec7 <- sort(sample(10^7, 10^3))

    ## 166,167,000 Combinations
    system.time(RcppAlgos::comboGeneral(tVec7, 3))
    system.time(arrangements::combinations(x = tVec7, k = 3))

    ## 124,251,000 Permuations
    system.time(RcppAlgos::permuteGeneral(tVec7[1:500], 3))
    system.time(arrangements::permutations(x = tVec7[1:500], k = 3))

    #' ## 7. Memory
    #'
    #' When executing `comboGeneral` as well as `arrangements::combinations`, the memory will jump almost 2 Gbs before calling `gc`. This seems about right as `#rows * #nols * bytesPerCell / 2^30 bytes = choose(1000,3) * 3 * 4 / 2^30 bytes = (166167000 * 3 * 4)/2^30 = 1.857 Gbs`). However, when executing `combn`, the memory behavior was eratic (e.g. sometimes it would use all 16 Gb of memory and other times it would only spike a couple of Gbs). When I tested this on the Windows set-up, it would often crash.
    #'
    #' We can confirm this using `Rprof` along with `summaryRporf`. Observe:

    Rprof("RcppAlgos.out", memory.profiling = TRUE)
    t1 <- RcppAlgos::comboGeneral(tVec7, 3)
    Rprof(NULL)
    head(summaryRprof("RcppAlgos.out", memory = "both")$by.total, n = 1)

    Rprof("arrangements.out", memory.profiling = TRUE)
    t3 <- arrangements::combinations(tVec7, 3)
    Rprof(NULL)
    head(summaryRprof("arrangements.out", memory = "both")$by.total, n = 1)

    #' With `RcppAlgos` & `arrangements`, `mem.total` registers just over `1900 Mb`.
    #'
    #' And here is the memory profile on a smaller vector.

    tVec7Prime <- tVec7[1:300]

    Rprof("combinat.out", memory.profiling = TRUE)
    t3 <- combinat::combn(tVec7Prime, 3)
    Rprof(NULL)
    head(summaryRprof("combinat.out", memory = "both")$by.total, n = 1)

    Rprof("utils.out", memory.profiling = TRUE)
    t4 <- utils::combn(tVec7Prime, 3)
    Rprof(NULL)
    head(summaryRprof("utils.out", memory = "both")$by.total, n = 1)

    Rprof("gtools.out", memory.profiling = TRUE)
    t5 <- gtools::combinations(300, 3, tVec7Prime)
    Rprof(NULL)
    head(summaryRprof("gtools.out", memory = "both")$by.total, n = 1)

    #' Interestingly, `utils::combn` and `combinat::combn` use different amounts of memory and take different amounts of time to execute.  This does not hold up with smaller vectors:

    microbenchmark(combinat::combn(2:13, 6), utils::combn(2:13, 6))

    #' And with `gtools` the total memory used is a little over 3x as much as `utils`. It should be noted that for these 3 packages, I obtained different results every-time I ran them (e.g. for `combinat::combn` sometimes I would get 9000 Mb and then I would get 13000 Mb).
    #'
    #' Still, none can match `RcppAlgos` **OR** `arrangements`. Both only use 51 Mb when ran on the example above.
    #'
    #' benchmark script: https://github.com/jwood000/RcppAlgos/blob/main/scripts/SO_Comb_Perm_in_R.R
    #'
    #' <sub>\*: An homage to _A Walk through Combinatorics_ by Miklós Bóna </sub>
}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
