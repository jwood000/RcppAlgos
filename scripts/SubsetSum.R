reprex::reprex({
    #'
    #' This document covers the topic of solving problems related to the [subset sum problem](<https://en.wikipedia.org/wiki/Subset_sum_problem>) with `RcppAlgos`. We have already covered integer partitions, which is a special case of the subset sum problem, in [Constraints and Integer Partitions](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#integer-partitions>) and it is highly encouraged to read that vignette first.
    #'
    #' ***
    #'
    #' ## Building on Integer Partitions
    #'
    #' The [integer partition problem](<https://en.wikipedia.org/wiki/Partition_\(number_theory\)>) presents the question _"how can we write n as a sum of positive integers?"_ There are well-known algorithms for enumerating _all_ partitions of an integer _n_. We even have algorithms for generating partitions of a specific length or with distinct parts only. But how do we enumerate partitions of _n_ with a specific set of numbers? What about enumerating partitions of a specific length _m_ of _n_ given a specific set of numbers?
    #'
    #' For example, using only the numbers `3:18`, find all partitions of `50` of length `5`.
    #'
    #' With `RcppAlgos`, this is easily achieved. We simply use the same template as we did in [Constraints and Integer Partitions](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#integer-partitions>). Observe (We continue to use the `ht` function defined in the [Combination and Permutation Basics](<https://jwood000.github.io/RcppAlgos/articles/GeneralCombinatorics.html>) vignette):
    #'

    library(RcppAlgos)
    options(width = 90)

    ht <- function(d, m = 5, n = m) {
        ## print the head and tail together
        cat("head -->\n")
        print(head(d, m))
        cat("--------\n")
        cat("tail -->\n")
        print(tail(d, n))
    }

    ## Each element can only occur once
    ht(partitionsGeneral(3:18, 5, target = 50))

    ## What about allowing repetition?
    ht(partitionsGeneral(3:18, 5, TRUE, target = 50))

    ## Even works on multisets
    ht(partitionsGeneral(3:18, 5, freqs = rep(1:4, 4), target = 50))

    #'
    #' In fact, these optimized algorithms can be applied when the vector passed has the quality that if you were to sort them, the difference of each element with it's neighbor is constant (E.g. `c(121, 126, 131, 136, ..., 221)`).
    #'

    even_time <- system.time({
        genParts <- partitionsGeneral(seq(121, 221, 5), 13,
                                      TRUE, target = 2613)
    })
    even_time

    ht(genParts)

    prettyNum(comboCount(seq(121, 221, 5), 13, TRUE), big.mark = ",")

    system.time(genMultiParts <- partitionsGeneral(seq(121, 221, 5), 13,
                                                   freqs = rep(1:7, 3),
                                                   targe = 2613))

    ht(genMultiParts)

    prettyNum(comboCount(seq(121, 221, 5), 13, freqs = rep(1:7, 3)), big.mark = ",")

    #' ### Working with Negative Numbers
    #'
    #' Generally, integer partition algorithms are restricted to positive integers. However, with the generalized partition algorithms in `RcppAlgos`, we can make light work of vectors with negative numbers (again, the sorted vector has to have the property that the difference of each element with it's neighbor is constant).
    #'

    system.time({
        genDistParts <- partitionsGeneral(seq(-173L, 204L, 13L),
                                          11, target = -460)
    })

    all(rowSums(genDistParts) == -460L)

    ht(genDistParts)

    #'
    #' ## Partitions with no Restrictions
    #'
    #' With the examples illustrated above, we had the restriction that the sorted input vector had to have the property that the difference of each element with it's neighbor is constant. If this requirement is broken, it only means that we cannot use a particular algorithm and we must fall back to a more general algorithm. _Fret not!!_ These general algorithms are extremely efficient and very flexible. We can use them with random input vectors, random targets, as well as over ranges.
    #'
    #' Let us revisit the example above but with a slight variation that breaks the requirement.
    #'

    inpVec <- c(116, seq(126, 221, 5))

    ## Non-constant difference... The specialized algo can't be used
    diff(inpVec)

    uneven_time <- system.time({
        genParts2 <- partitionsGeneral(inpVec, 13, TRUE, target = 2613)
    })
    uneven_time    ## out of a possible 573 million in under a second

    ht(genParts2)

    #'
    #' Although the above was about `r round(uneven_time[["elapsed"]] / even_time[["elapsed"]])` times slower than the first example dealing with 573 million combinations (`r uneven_time[["elapsed"]]` milliseconds vs. `r even_time[["elapsed"]]` milliseconds), we are still dealing in milliseconds!!! For reference, version `2.3.4` takes about 18 seconds to find all 118,560 solutions.
    #'
    #' Here are some more exotic examples demonstrating the power of these algorithms.
    #'

    set.seed(42)
    mySamp <- sample(-100:100, 50)

    sort(mySamp)

    system.time(exotic <- partitionsGeneral(mySamp, 8, freqs = rep(1:5, 10),
                                            target = 496))

    dim(exotic)

    ## Over 1 billion total combinations
    prettyNum(comboCount(mySamp, 8, freqs = rep(1:5, 10)), big.mark = ",")

    ## Only getting a few (a thousand in this case) is much faster
    system.time(partitionsGeneral(mySamp, 8, freqs = rep(1:5, 10),
                                  target = 496, upper = 1e3))

    #'
    #' The function `permuteGeneral` benefits from these optimized algorithms as well. However, just as we discussed in [Output Order with `permuteGeneral`](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#output-order-with-permutegeneral>), the output will not be in lexicographical order.
    #'
    #' ## Taming Floating Point Numbers
    #'
    #' Oftentimes when working with numerical vectors, it can be hard to find combinations that sum to a particular number because of floating point errors (See [Using `tolerance`](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html#using-tolerance>)). In practice, we may not need an exact match and a close approximation will suffice. For example, let's say we have a football team of 100 players (including practice squad) and we are interested in a trade involving 6 players and a total salary of 20 million dollars. We may not be able to find 6 players whose sum of salaries is exactly 20 million, but we can find many 6 player combinations whose sum of salaries is within a tolerance of 20 million.
    #'

    set.seed(22213)
    football_player_salaries <- 2e7 * rbeta(100, 2, 25)

    summary(football_player_salaries)

    ## Over 1 billion combinations...
    ## An exhaustive search will not be feasible
    prettyNum(comboCount(football_player_salaries, 6), big.mark = ",")

    system.time(exactly20 <- partitionsGeneral(football_player_salaries, 6,
                                               target = 2e7, tolerance = 0))

    ## No results that equal exactly 2e7
    dim(exactly20)

    #'
    #' What if we increase the tolerance to `$`1000 (Honestly... what's `$`1000 when we are talking about 20 million dollars)? Our intent is to explore these options, so we take advantage of the `upper` argument in anticipation that we obtain many results that meet the criteria. If we obtain the upper bound, we decrease the tolerance (if needed) and repeat.
    #'

    ## N.B. This is much more efficient. Also, we set keepResults
    ## to TRUE so we can see the total sum of salaries.
    system.time(almost20 <- comboGeneral(football_player_salaries, 6,
                                         constraintFun = "sum", comparisonFun = "==",
                                         limitConstraints = 2e7, tolerance = 1000,
                                         upper = 1000, keepResults = TRUE))

    dim(almost20)

    ht(almost20)

    ## decreasing the tolerance to $10 further we obtain 158 results
    system.time(superClose20 <- comboGeneral(football_player_salaries, 6,
                                             constraintFun = "sum", comparisonFun = "==",
                                             limitConstraints = 2e7, tolerance = 10,
                                             upper = 1000, keepResults = TRUE))

    ht(superClose20)

    #'
    #' ## `prod` and `mean`
    #'
    #' These optimized algorithms are also employed when `constraintFun` is `"prod"` or `"mean"`.
    #'

    getAllThenFilter <- function(n, m, lim) {
        t <- comboGeneral(n, m, constraintFun = "prod")
        t[t[, m + 1] == lim, -(m+1)]
    }

    library(microbenchmark)
    microbenchmark(optimized = comboGeneral(25, 10, constraintFun = "prod",
                                            comparisonFun = "==",
                                            limitConstraints = 1037836800),
                   brute = getAllThenFilter(25, 10, 1037836800), times = 20,
                   unit = "relative", check = "equal")

    ## What about cases when brute force isn't feasible
    set.seed(101)
    v <- runif(1000, 1, 2)

    prettyNum(comboCount(v, 100), big.mark = ",")

    system.time(prodAlmost100 <- comboGeneral(v, 100, constraintFun = "prod",
                                              comparisonFun = "==",
                                              limitConstraints = 100,
                                              tolerance = 0.0001, upper = 20))

    dim(prodAlmost100)

    apply(prodAlmost100, 1, prod)

    ## Showcasing mean
    system.time(meanAlmost1.5 <- comboGeneral(v, 100, constraintFun = "mean",
                                              comparisonFun = "==",
                                              limitConstraints = 1.5,
                                              tolerance = 0.0001, upper = 20))

    dim(meanAlmost1.5)

    rowMeans(meanAlmost1.5)

    #'
    #' ## Using Iterators
    #'
    #' As of version `2.5.0` all of the above cases can be attacked with iterators (See [Combinatorial Iterators in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/CombinatoricsIterators.html#iterating-over-constrained-combinations/permutations>)). As mentioned in the suggested reading, iterators are very flexible and just as efficient as their "general" counterparts. They have the added benefit of allowing one to save the current state, allowing one to generate _n_ results at a time.
    #'
    #' Below are a few demonstrations using some of the examples in earlier sections.
    #'

    ## The football salary example
    salary <- partitionsIter(football_player_salaries, 6,
                             target = 2e7, tolerance = 1000)

    ## Or use comboIter:
    ##
    ## comboIter(football_player_salaries, 6, constraintFun = "sum",
    ##           comparisonFun = "==", limitConstraints = 2e7,
    ##           tolerance = 1000, upper = 1000, keepResults = TRUE))

    system.time(almost20withIter <- salary@nextNIter(1e3))

    ## almost20 was generated above with comboGeneral
    all.equal(almost20[, 1:6], almost20withIter)

    ## With iterators we can easily continue iterating. With the general
    ## functions if we wanted the next 1000 results, we would have to
    ## generate the first 1000 along with the next 1000
    system.time(nextAlmost20withIter <- salary@nextNIter(1e3))

    ht(nextAlmost20withIter)

    salary@summary()


    ## The prodAlmost100 example
    prodIter <- comboIter(v, 100,
                          constraintFun = "prod", comparisonFun = "==",
                          limitConstraints = 100, tolerance = 0.0001)

    system.time(prodAlmost100WithIter <- prodIter@nextNIter(20))

    all.equal(prodAlmost100, prodAlmost100WithIter)

    ## Again, with iterators, we can continue iterating from
    ## where we left off
    system.time(nextAlmost100WithIter <- prodIter@nextNIter(20))

    dim(nextAlmost100WithIter)

    ## Use @ or $ to access methods. If one needs to access these methods
    ## often (e.g. nextIter inside a loop), it is recommended to use the
    ## @ accessor as it is much more efficient.
    prodIter$summary()

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")