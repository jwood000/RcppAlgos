reprex::reprex({
    #'
    #' This document covers topics in generating random samples of combinations, permutations, partitions, compositions and partition of groups. It is encouraged to read [General Combinatorics](<https://jwood000.github.io/RcppAlgos/articles/GeneralCombinatorics.html>) first.
    #'
    #' ***
    #'
    #' # Sampling
    #'
    #' At the heart of sampling is the ability to efficiently generate the _n<sup>th</sup>_ [lexicographical](<https://en.wikipedia.org/wiki/Lexicographic_order>) result. The algorithms in `RcppAlgos` are flexible and optimized, allowing for tackling this task with ease.
    #'
    #' ## Base R
    #'
    #' To illustrate this in base R, let us consider getting 5 random combinations of the vector `1:20` of length 10. How should we proceed?
    #'
    #' A naive approach would be to generate all of the combinations using `combn` and then call `sample`:
    #'

    options(width = 90)
    naive <- function(v, m, n, s) {
        allCombs <- combn(v, m)
        set.seed(s)
        allCombs[, sample(ncol(allCombs), n)]
    }

    fiveRndCombs <- naive(20, 10, 5, 42)
    t(fiveRndCombs)

    #'
    #' This is okay for this small example (there are only `choose(20, 10) = 184756` results), however what if we wanted to find one hundred thousand random combinations from the vector `1:100` of length 20? Clearly, the approach above will not be feasible as there are far too many results to generate (`choose(100, 20) = 5.359834e+20`). Furthermore, there are internal limitations on `sample`. If we try to pass `choose(100, 20)`, we will get an error:
    #'

    sample(choose(100, 20), 5)

    #'
    #'
    #' We could also try calling `sample(100, 20)` a bunch of times and hope we don't get duplicate combinations. This is neither promising nor elegant.
    #'
    #' ***
    #'
    #' ## RcppAlgos Solutions
    #'
    #' `RcppAlgos` provides five functions: `comboSample`, `permuteSample`, `partitionsSample`, `compositionsSample`, and `comboGroupsSample` for seamlessly attacking these types of problems. All functions provide the following:
    #'
    #' * Easily generate random samples in parallel using the `nThreads` or the `Parallel` parameters.
    #' * You can pass a vector of specific indices via `sampleVec` or rely on the internal sampling functions. We call `sample` when the total number of results is small and for larger cases, the sampling is done in a very similar fashion to `urand.bigz` from the `gmp` package.
    #' * Consistent interface to their respective general functions (e.g. `partitionsGeneral`)
    #' * The `seed` parameter allows for generating reproducible samples.
    #' * If the gmp library is needed, the `seed` parameter must be set in order to have reproducible results (_E.g._ `set.seed()`) has no effect in these cases).
    #'
    #' ## `comboSample` and `permuteSample`
    #'
    #' Let's first look at the first problem above (i.e. getting 5 random combinations of the vector `1:20` of length 10):
    #'

    library(RcppAlgos)
    set.seed(42)
    comboSample(20, 10, n = 5)

    ## Use the seed argument directly to produce the same output
    comboSample(20, 10, n = 5, seed = 42)

    ## fiveRndCombs produced above
    identical(t(fiveRndCombs),
              comboSample(20, 10, n = 5, seed = 42))

    #'
    #'
    #' ### Samples of Results with Repetition
    #'
    #' Just like with `comboGeneral` and `permuteGeneral`, we can explore results with repetition.
    #'

    comboSample(10, 8, TRUE, n = 3, seed = 84)

    permuteSample(10, 8, TRUE, n = 3)

    comboSample(10, 12, freqs = 1:10, n = 3)

    permuteSample(10, 12, freqs = 1:10, n = 3, seed = 123)

    #'
    #' ### Specific Results with `sampleVec`
    #'
    #' We can also utilize `sampleVec` to generate specific results.
    #'

    ## E.g. the below generates the 1st, 5th, 25th, 125th, and
    ## 625th lexicographical combinations
    comboSample(10, 8, TRUE, sampleVec = c(1, 5, 25, 125, 625))

    ## Is the same as:
    comboGeneral(10, 8, TRUE)[5^(0:4), ]

    #'
    #' ### Using `namedSample`
    #'
    #' Have you ever wondered which lexicographical combinations/permutations are returned when sampling? No worries, simply set `namedSample = TRUE`:
    #'

    testInd <- permuteSample(30, 10, n = 3, seed = 100, namedSample = TRUE)
    testInd

    ## Same output as above
    permuteSample(30, 10, sampleVec = row.names(testInd))

    #'
    #' ### Parallel Computing and GMP Support
    #'
    #' Just like the `General` counterparts, the sampling functions utilize GMP to allow for exploration of combinations/permutations of large vectors where the total number of results is enormous. They also offer parallel options using `Parallel` or `nThreads`.
    #'

    ## Uses min(stdThreadMax() - 1, 5) threads (in this case)
    permuteSample(500, 10, TRUE, n = 5, seed = 123, Parallel = TRUE)

    permuteSample(factor(state.abb), 15, n = 3, seed = 50, nThreads = 3)

    permuteCount(factor(state.abb), 15)

    #'
    #' ### Efficiency
    #'
    #' The algorithms are incredibly efficient and offer tremendous gains over the naive approach above:
    #'

    ## the function "naive" is defined above
    system.time(naive(25, 10, 5, 15))

    system.time(comboSample(25, 10, n = 5, seed = 15))

    #'
    #' Even when dealing with extremely large numbers, these algorithms are very fast. And using the parallel options have even greater effects than we saw with the general counterparts (typically around ~2-3 times faster with the general functions, whereas with the last example below with sampling we see a nearly 5x improvement).
    #'

    ## Lightning fast even with examples involving many results
    system.time(comboSample(2500, 100, n = 5, seed = 15))

    ## The total number of combinations has ~180 digits
    gmp::log10.bigz(comboCount(2500, 100))

    ## Still fast with larger samples
    system.time(comboSample(2500, 100, n = 1e4, seed = 157))

    ## Using Parallel/nThreads in these cases has an even greater effect
    system.time(comboSample(2500, 100, n = 1e4, seed = 157, nThreads = 8))

    #'
    #' ### User Defined Functions
    #'
    #' Again, just as with the general functions, you can pass a custom function to `{combo|permute}Sample` using the `FUN` argument.
    #'

    permuteSample(5000, 1000, n = 3, seed = 101, FUN = sd)

    ## Example using complex numbers
    myCplx <- as.complex(1:100 + rep(c(-1, 1), 50) * 1i)

    permuteSample(myCplx, 10, freqs = rep(1:5, 20),
                  n = 3, seed = 101, FUN = function(x) {
                      sqrt(sum(x))
                  })

    #'
    #' ## `partitionsSample`
    #'
    #' The `partitionsSample` function allows one to draw a random sample of partitions of a number. Many of the features present in `comboSample` and `permuteSample` are available in `partitionsSample`.
    #'

    ## Use the seed parameter to obtain reproducible results
    partitionsSample(100, 8, TRUE, n = 3, seed = 42)

    ## Used namedSample to obtain the lexicographical indices
    partitionsSample(100, 8, TRUE, n = 3, seed = 42, namedSample = TRUE)

    ## Use sampleVec to obtain specific results
    partitionsSample(100, 8, TRUE, sampleVec = c(61413, 54425, 623844))

    partitionsCount(2500, 10)

    ## Algorithms are very efficient
    system.time(serial <- partitionsSample(2500, 10, n = 1e3,
                                           seed = 8128))

    ## Use nThreads for greater efficiency
    system.time(multi <- partitionsSample(2500, 10, n = 1e3,
                                          seed = 8128, nThreads = 8))

    identical(multi, serial)

    ## Even works with non-standard setup
    partitionsSample(17 + (1:10) * 3, 10, TRUE,
                     target = 320, n = 3, seed = 111)

    #'
    #' There are sampling algorithms available for most partition cases, but some cases are not covered. For example, with standard multisets, we are currently unable to _efficiently_ generate the _n<sup>th</sup>_ lexicographical result. Another example is when the source vector is not uniform (_e.g._ when the distance between each element is irregular).
    #'
    #' Observe the following:
    #'

    ## No sampling algorithm available when the source vector is not uniform
    partitionsSample(c(1, 4, 6, 7, 10, seq(11, 100, 7)), 10, n = 1, target = 340)

    ## As stated above, the standard multiset case doesn't work either
    partitionsSample(0:50, 6, freqs = rep(1:3, 17), n = 2)

    ## If we use freqs to indicate that zeros can repeat,
    ## then we can obtain random samples
    partitionsSample(0:50, 6, freqs = c(50, rep(1, 50)), n = 3, seed = 222)

    ## Even works when the vector is restricted in regards to the target
    partitionsSample(0:50, 6, freqs = c(50, rep(1, 50)),
                     n = 3, seed = 222, target = 100)
    #'
    #'
    #' There is ongoing research in this area and our goal is to eventually be able to cover the standard multiset case.
    #'
    #' ## `compositionsSample`
    #'
    #' The `compositionsSample` function allows one to draw a random sample of compositions of a number. Many of the features present in `comboSample` and `permuteSample` are available in `compositionsSample`.
    #'

    ## Use the seed parameter to obtain reproducible results
    compositionsSample(100, 8, TRUE, n = 3, seed = 42)

    ## Used namedSample to obtain the lexicographical indices
    compositionsSample(100, 8, TRUE, n = 3, seed = 42, namedSample = TRUE)

    ## Use sampleVec to obtain specific results
    compositionsSample(100, 8, TRUE,
                       sampleVec = c(4024715585, 2756281572, 4873365553))

    compositionsCount(2500, 10, TRUE)

    ## Algorithms are very efficient...
    ## The below retrieves 10,000 compositions in under a second
    system.time(serial <- compositionsSample(2500, 10, TRUE,
                                             n = 1e4, seed = 8128))

    ## Use nThreads for greater efficiency
    system.time(multi <- compositionsSample(2500, 10, TRUE, n = 1e4,
                                            seed = 8128, nThreads = 8))

    identical(multi, serial)

    ## Sample weak compositions
    compositionsSample(0:100, 8, repetition = TRUE, weak = TRUE,
                       seed = 245659, n = 3, namedSample = TRUE)
    #'
    #'
    #' Currently, there are only sampling algorithms for most cases of compositions with repetition. There is ongoing work to expand these algorithms in the future.
    #'
    #' ## Sampling Partitions of Groups of Equal Size with `comboGroupsSample`
    #'
    #' Just as we can generate random samples of combinations and permutations, we are also able to generate random samples of partitions of groups of equal size. There are many problems that present in this manner. Below, we examine one involving playing cards.
    #'
    #' Let's say we have 4 players and each player is to have 3 cards a piece. Given that the deck is shuffled, the dealer then distrubutes 12 cards.
    #'
    #' > What possible hands can each player have?
    #'
    #' See [Creating A Deck Of Cards In R Without Using While And Double For Loop
    #' ](https://stackoverflow.com/a/36903806/4408538) (Credit to @MichaelChirico)
    #'

    cards <- c(2:10, "J", "Q", "K", "A")
    suits <- c("♠", "♥", "♦", "♣")
    deck <- paste0(rep(cards, length(suits)),  #card values
                   rep(suits, each = length(cards))) #suits

    set.seed(1738)
    shuffled <- factor(deck[sample(52)], levels = deck)

    ## Here are 3 possibilities
    comboGroupsSample(shuffled[1:12], numGroups = 4, n = 2, seed = 13)


    comboGroupsSample(shuffled[1:12], numGroups = 4, retType = "3Darray",
                      n = 2, seed = 13, namedSample = TRUE)

    #'
    #' # Ranking
    #'
    #' Ranking is the complement of sampling. That is, given a combination (or permutation/partition), determine which lexicographical result it is. As an example, consider all of the combinations of 5 choose 3:
    #'

    comboGeneral(5, 3)

    #'
    #' We can see that the rank of the combination: `c(2, 3, 4)` is 7. That is, `c(2, 3, 4)` is the 7<sup>th</sup> combination of 5 choose 3.
    #'
    #' ## Base R
    #'
    #' Just as we saw before, we could easily produce a brute force approach that would work well with small cases, but would become unmanagemable very quickly. For example:
    #'

    naive_rank <- function(v, m, comb) {
        comb <- as.integer(comb)
        which(apply(combn(v, m), 2, function(x) identical(x, comb)))
    }

    naive_rank(5, 3, 2:4)

    ## Larger example
    comb = comboSample(25, 12, sampleVec = 2e6)[1, ]

    system.time(print(naive_rank(25, 12, comb)))
    #'
    #'
    #' ## RcppAlgos Solutions
    #'
    #' Similar to the sampling problem, `RcppAlgos` provides four functions: `comboRank`, `permuteRank`, `partitionsRank`, and `compositionsRank` (currently there is not a ranking function for `comboGroups`). These functions are very similar to their sampling counterparts.
    #'
    #' For both problems presented above, here is how you would attack them with `comboRank`:
    #'

    comboRank(2:4, v = 5)

    ## Since order doesn't matter with combinations, 4:2 should return 7 as well
    comboRank(4:2, v = 5)

    ## comb was provided above
    system.time(print(comboRank(comb, v = 25)))
    #'
    #'
    #' All that is needed is the original vector that was used to produce the results and whether or not repetition is used via the `repetition` or `freqs` arguments. The width is determined automatically by the input.
    #'
    #' ## Rank Multiple Inputs
    #'
    #' A neat feature of the ranking functions is the ability to rank multiple inputs at once. We can either pass a single vector, multiple vectors, and we can even pass matrices.
    #'
    #' ## `comboRank`
    #'

    combs = comboSample(50, 8, n = 10, seed = 123, namedSample = TRUE)
    combs

    comboRank(combs, v = 50)
    #'
    #'
    #' ## `permuteRank`
    #'

    perms_len_5 = permuteSample(100, 5, freqs = rep(1:5, 20),
                                n = 3, seed = 987, namedSample = TRUE)
    perms_len_5

    perms_len_8 = permuteSample(100, 8, freqs = rep(1:5, 20),
                                n = 3, seed = 123, namedSample = TRUE)
    perms_len_8

    ## Note you can name the inputs
    permuteRank(p5 = perms_len_5, p8 = perms_len_8,
                v = 100, freqs = rep(1:5, 20))
    #'
    #'
    #' ## `partitionsRank`
    #'

    parts = partitionsSample(50, 8, target = 100, repetition = TRUE,
                             n = 3, seed = 42, namedSample = TRUE)
    parts

    partitionsRank(parts, v = 50, target = 100, repetition = TRUE)
    #'
    #'
    #' ## `compositionsRank`
    #'

    comps = compositionsSample(50, 8, repetition = TRUE,
                               n = 3, seed = 42, namedSample = TRUE)
    comps

    compositionsRank(comps, v = 50, repetition = TRUE)

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")