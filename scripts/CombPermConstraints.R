reprex::reprex({
    #'
    #' This document covers the topic of finding combinations or permutations that meet a specific set of criteria. For example, retrieving all combinations of a vector that have a product between two bounds.
    #'
    #' ***
    #'
    #' ## Constraint Functions
    #'
    #' There are 5 compiled constraint functions that can be utilized efficiently to test a given result.
    #'
    #'   1. sum
    #'   2. prod
    #'   3. mean
    #'   4. max
    #'   5. min
    #'
    #' They are passed as strings to the `constraintFun` parameter. When these are employed without any other parameters being set, an additional column is added that represents the result of applying the given function to that combination/permutation. You can also set `keepResults = TRUE` (more on this later).
    #'

    library(RcppAlgos)
    options(width = 90)

    packageVersion("RcppAlgos")

    cat(paste(capture.output(sessionInfo())[1:3], collapse = "\n"))

    ## base R using combn and FUN
    combnSum = combn(20, 10, sum)
    algosSum = comboGeneral(20, 10, constraintFun = "sum")

    ## Notice the additional column (i.e. the 11th column)
    head(algosSum)

    identical(as.integer(combnSum), algosSum[,11])

    ## Using parallel
    paralSum = comboGeneral(20, 10, constraintFun = "sum", Parallel = TRUE)
    identical(paralSum, algosSum)

    library(microbenchmark)
    microbenchmark(serial = comboGeneral(20, 10, constraintFun = "sum"),
                 parallel = comboGeneral(20, 10, constraintFun = "sum", Parallel = TRUE),
                 combnSum = combn(20, 10, sum), unit = "relative")

    #'
    #' ### Faster than `rowSums` and `rowMeans`
    #'
    #' Finding row sums or row means is even faster than simply applying the highly efficient `rowSums`/`rowMeans` _after_ the combinations have already been generated:
    #'

    ## Pre-generate combinations
    combs = comboGeneral(25, 10)

    ## Testing rowSums alone against generating combinations as well as summing
    microbenchmark(serial = comboGeneral(25, 10, constraintFun = "sum"),
                 parallel = comboGeneral(25, 10, constraintFun = "sum", Parallel = TRUE),
                  rowsums = rowSums(combs), unit = "relative")

    all.equal(rowSums(combs),
              comboGeneral(25, 10,
                           constraintFun = "sum",
                           Parallel = TRUE)[,11])

    ## Testing rowMeans alone against generating combinations as well as obtain row means
    microbenchmark(serial = comboGeneral(25, 10, constraintFun = "mean"),
                 parallel = comboGeneral(25, 10, constraintFun = "mean", Parallel = TRUE),
                 rowmeans = rowMeans(combs), unit = "relative")

    all.equal(rowMeans(combs),
              comboGeneral(25, 10,
                           constraintFun = "mean",
                           Parallel = TRUE)[,11])

    #'
    #' In both cases above, `RcppAlgos` is doing double the work nearly twice as fast!!!
    #'
    #' ## Comparison Operators and `limitConstraints`
    #'
    #' The standard 5 comparison operators (i.e. `"<"`, `">"`, `"<="`, `">="`, & `"=="`) can be used in a variety of ways. In order for them to have any effect, they must be used in conjunction with `constraintFun` as well as `limitConstraints`. The latter is the value(s) that will be used for comparison. It can be passed as a single value or a vector of two numerical values. This is useful when one wants to find results that are between (or outside) of a given range.
    #'
    #' ### One Comparison Operator
    #'
    #' First we will look at cases with only one comparison and one value for the `limitConstraint`.
    #'

    ## Generate some random data. N.B. Using R >= 4.0.0
    set.seed(101)
    myNums = sample(500, 20)

    myNums

    ## Find all 5-tuples combinations without repetition of myNums
    ## (defined above) such that the sum is equal to 1176.
    p1 = comboGeneral(v = myNums, m = 5,
                      constraintFun = "sum",
                      comparisonFun = "==",
                      limitConstraints = 1176)

    tail(p1)


    ## Authenticate with brute force
    allCombs = comboGeneral(sort(myNums), 5)
    identical(p1, allCombs[which(rowSums(allCombs) == 1176), ])


    ## How about finding combinations with repetition
    ## whose mean is less than or equal to 150.
    p2 = comboGeneral(v = myNums, m = 5, TRUE,
                      constraintFun = "mean",
                      comparisonFun = "<=",
                      limitConstraints = 150)

    ## Again, we authenticate with brute force
    allCombs = comboGeneral(sort(myNums), 5, TRUE)
    identical(p2, allCombs[which(rowMeans(allCombs) <= 150), ])

    ## N.B.
    class(p2[1, ])

    class(allCombs[1, ])

    ## When mean is employed or it can be determined that integral
    ## values will not suffice for the comparison, we fall back to
    ## numeric types, thus all.equal should return TRUE
    all.equal(p2, allCombs[which(rowMeans(allCombs) <= 150), ])

    #'
    #' ### Two Comparison Operators
    #'
    #' Sometimes, we need to generate combinations/permutations such that when we apply a constraint function, the results are between (or outside) a given range. There is a natural two step process when finding results outside a range, however for finding results between a range, this two step approach could become computationally demanding. The underlying algorithms in `RcppAlgos` are optimized for both cases and avoids adding results that will eventually be removed.
    #'
    #' Using two comparisons is easy. The first comparison operator is applied to the first limit and the second operator is applied to the second limit.
    #'
    #' Note that in the examples below, we have `keepResults = TRUE`. This means an additional column will be added to the output that is the result of applying `constraintFun` to that particular combination.
    #'

    ## Get combinations such that the product is
    ## strictly between 3600 and 4000
    comboGeneral(5, 7, TRUE, constraintFun = "prod",
                 comparisonFun = c(">","<"),          ## Find results > 3600 and < 4000
                 limitConstraints = c(3600, 4000),
                 keepResults = TRUE)

    # ## The above is the same as doing the following:
    # comboGeneral(5, 7, TRUE, constraintFun = "prod",
    #              comparisonFun = c("<",">"),          ## Note that the comparison vector
    #              limitConstraints = c(4000, 3600),    ## and the limits have flipped
    #              keepResults = TRUE)


    ## What about finding combinations outside a range
    outside = comboGeneral(5, 7, TRUE, constraintFun = "prod",
                           comparisonFun = c("<=",">="),
                           limitConstraints = c(3600, 4000),
                           keepResults = TRUE)

    all(apply(outside[, -8], 1, prod) <= 3600
         | apply(outside[, -8], 1, prod) >= 4000)

    dim(outside)

    ## Note that we obtained 5 results when searching "between"
    ## 3600 and 4000. Thus we have: 325 + 5 = 330
    comboCount(5, 7, T)

    #'
    #' ### Using `tolerance`
    #'
    #' When the underlying type is `numeric`, [round-off errors](<https://en.wikipedia.org/wiki/Round-off_error>) can occur. As stated in [floating-point error mitigation](<https://en.wikipedia.org/wiki/Floating-point_error_mitigation>):
    #'
    #' > _**"By definition, floating-point error cannot be eliminated, and, at best, can only be managed."**_
    #'
    #' Here is a great stackoverflow post that further illuminates this tricky topic:
    #'
    #'   - [What is the correct/standard way to check if difference is smaller than machine precision?](<https://stackoverflow.com/q/59229545/4408538>)
    #'
    #' For these reasons, the argument `tolerance` can be utilized to refine a given constraint. It is added to the upper limit and subtracted from the lower limit. The default value is `sqrt(.Machine$double.eps) ~= 0.00000001490116`.
    #'
    #' This default value is good and bad.
    #'
    #' For the good side:
    #'

    dim(comboGeneral(seq(0, 0.5, 0.05), 6, TRUE,
                     constraintFun = "sum",
                     comparisonFun = "==",
                     limitConstraints = 1))

    ## Confirm with integers and brute force
    allCbs = comboGeneral(seq(0L, 50L, 5L), 6, TRUE, constraintFun = "sum")

    sum(allCbs[, 7] == 100L)

    #'
    #' If we had a tolerance of zero, we would have obtained an incorrect result:
    #'

    ## We miss 31 combinations that add up to 1
    dim(comboGeneral(seq(0, 0.5, 0.05), 6, TRUE,
                     constraintFun = "sum",
                     comparisonFun = "==",
                     limitConstraints = 1, tolerance = 0))

    #'
    #' And now for a less desirable result. The example below appears to give incorrect results. That is, we shouldn't return any combination with a mean of 4.1 or 5.1.
    #'

    comboGeneral(c(2.1, 3.1, 5.1, 7.1), 3, T,
                 constraintFun = "mean", comparisonFun = c("<", ">"),
                 limitConstraints = c(5.1, 4.1), keepResults = TRUE)

    #'
    #' In the above example, the range that is actually tested against is `c(4.0999999950329462, 5.1000000049670531)`.
    #'
    #' If you want to be absolutely sure you are getting the correct results, one must rely on integers as simple changes in arithmetic can throw off precision in floating point operations.
    #'

    comboGeneral(c(21, 31, 51, 71), 3, T,
                 constraintFun = "mean", comparisonFun = c("<", ">"),
                 limitConstraints = c(51, 41), keepResults = TRUE) / 10

    #'
    #' ### Output Order with `permuteGeneral`
    #'
    #' Typically, when we call `permuteGeneral`, the output is in lexicographical order, however when we apply a constraint, the underlying algorithm checks against combinations only, as this is more efficient. If a particular combination meets a constraint, then all permutations of that vector also meet that constraint, so there is no need to check them. For this reason, the output isn't in order. Observe:
    #'

    permuteGeneral(c(2, 3, 5, 7), 3, freqs = rep(2, 4),
                   constraintFun = "mean", comparisonFun = c(">", "<"),
                   limitConstraints = c(4, 5), keepResults = TRUE, tolerance = 0)

    #'
    #' As you can see, the _2<sup>nd</sup>_ through the _6<sup>th</sup>_ entries are simply permutations of the _1<sup>st</sup>_ entry. Similarly, entries _8_ and _9_ are permutations of the _7<sup>th</sup>_ and entries _11_ and _12_ are permutations of the _10<sup>th</sup>_.
    #'
    #' ## Integer Partitions
    #'
    #' Specialized algorithms are employed when it can be determined that we are looking for [integer partitions](<https://en.wikipedia.org/wiki/Partition_(number_theory)>).
    #'
    #' As of version `2.5.0`, we now have added `partitionsGeneral` which is similar to `comboGeneral` with `constraintFun = "sum"` and `comparisonFun = "=="`. Instead of using the very general `limitConstraints` parameter, we use `target` with a default of `max(v)` as it seems more fitting for partitions.
    #'
    #' ### Case 1: All Integer Partitions of _N_
    #'
    #' We need `v = 0:N`, `repetition = TRUE`. When we leave `m = NULL`, `m` is internally set to the length of the longest non-zero combination (this is true for all cases below).
    #'

    partitionsGeneral(0:5, repetition = TRUE)

    ## Note that we could also use comboGeneral:
    ## comboGeneral(0:5, repetition = TRUE,
    ##              constraintFun = "sum",
    ##              comparisonFun = "==", limitConstraints = 5)
    ##
    ## The same goes for any of the examples below

    #'
    #' ### Case 2: Integer Partitions of _N_ of Length _m_
    #'
    #' Everything is the same as above except for explicitly setting the desired length and deciding whether to include zero or not.
    #'

    ## Including zero
    partitionsGeneral(0:5, 3, repetition = TRUE)

    ## Zero not included
    partitionsGeneral(5, 3, repetition = TRUE)

    #'
    #' ### Case 3: Integer Partitions of _N_ into Distinct Parts
    #'
    #' Same as `Case 1 & 2` except now we have `repetition = FALSE`.
    #'

    partitionsGeneral(0:10)

    ## Zero not included and restrict the length
    partitionsGeneral(10, 3)

    ## Include zero and restrict the length
    partitionsGeneral(0:10, 3)

    ## partitions of 10 into distinct parts of every length
    lapply(1:4, function(x) {
        partitionsGeneral(10, x)
    })

    #'
    #' #### Using `freqs` to Refine Length
    #'
    #' We can utilize the `freqs` argument to obtain more distinct partitions by allowing for repeated zeros. The super optimized algorithm will only be carried out if zero is included and the number of repetitions for every number except zero is one.
    #'
    #' For example, given `v = 0:N` and `J >= 1`, if `freqs = c(J, rep(1, N))`, then the super optimized algorithm will be used, however if `freqs = c(J, 2, rep(1, N - 1))`, the general algorithm will be used. It should be noted that the general algorithms are still highly optimized so one should not fear using it.
    #'
    #' A pattern that is guaranteed to retrieve all distinct partitions of _N_ is to set `v = 0:N` and `freqs = c(N, rep(1, N))` (the extra zeros will be left off).
    #'

    ## Obtain all distinct partitions of 10
    partitionsGeneral(0:10, freqs = c(10, rep(1, 10)))    ## Same as c(3, rep(1, 10))

    #' #### Caveats Using `freqs`
    #'
    #' As noted in `Case 1`, if `m = NULL`, the length of the output will be determined by the longest non-zero combination that sums to _N_.
    #'

    ## m is NOT NULL and output has at most 2 zeros
    partitionsGeneral(0:10, 3, freqs = c(2, rep(1, 10)))

    ## m is NULL and output has at most 2 zeros
    partitionsGeneral(0:10, freqs = c(2, rep(1, 10)))

    #'
    #' ### Case 4: Integer Partitions of _N_ into Parts of Varying Multiplicity
    #'

    ## partitions of 12 into 4 parts where each part can
    ## be used a specific number of times (e.g. 2 or 3)
    partitionsGeneral(12, 4, freqs = rep(2:3, 6))

    #'
    #' ## Efficiency Generating Partitions
    #'
    #' Note, as of version `2.5.0`, one can generate partitions in parallel using the `nThreads` argument.
    #'

    ## partitions of 60
    partitionsCount(0:60, repetition = TRUE)

    ## Single threaded
    system.time(partitionsGeneral(0:60, repetition = TRUE))

    ## Using nThreads
    system.time(partitionsGeneral(0:60, repetition = TRUE, nThreads=4))


    ## partitions of 120 into distinct parts
    partitionsCount(0:120, freqs = c(120, rep(1, 120)))

    system.time(partitionsGeneral(0:120, freqs = c(120, rep(1, 120))))

    system.time(partitionsGeneral(0:120, freqs = c(120, rep(1, 120)), nThreads=4))


    ## partitions of 100 into parts of 15 with specific multiplicity
    partitionsCount(100, 15, freqs = rep(4:8, 20))

    ## Over 6 million in just over a second!
    system.time(partitionsGeneral(100, 15, freqs = rep(4:8, 20)))

    #'
    #' ## Integer Compositions
    #'
    #' [Compositions](<https://en.wikipedia.org/wiki/Composition_(combinatorics)>) are related to integer partitions, however order matters. With `RcppAlgos`, we generate standard compositions with `compositionsGeneral`. Currently, the composition algorithms are limited to a subset of cases of compositions with repetiion.
    #'
    #' The output with `compositionGeneral` will be in lexicographical order. When we set `weak = TRUE`, we will obtain **_weak compositions_**, which allow for zeros to be a part of the sequence (E.g. `c(0, 0, 5), c(0, 5, 0), c(5, 0, 0)` are weak compositions of 5). As the Wikipedia article points out, we can increase the number of zeros indefinitely when `weak = TRUE`.
    #'
    #' For more general cases, we can make use of `permuteGeneral`, keeping in mind that the output will not be in lexicographical order. Another consideration with `permuteGeneral` is that when we include zero, we will always obtain weak compositions.
    #'
    #' With that in mind, generating compositions with `RcppAlgos` is easy, flexible, and quite efficient.
    #'
    #' ### Case 5: All Compositions of _N_
    #'

    ## See Case 1
    compositionsGeneral(0:3, repetition = TRUE)

    ## Get weak compositions
    compositionsGeneral(0:3, repetition = TRUE, weak = TRUE)

    ## Get weak compositions with width > than target
    tail(compositionsGeneral(0:3, 10, repetition = TRUE, weak = TRUE))

    ## With permuteGeneral, we always get weak compositions, just
    ## not in lexicographical order
    permuteGeneral(0:3, repetition = TRUE,
                   constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 3)

    tail(permuteGeneral(0:3, 10, repetition = TRUE,
                        constraintFun = "sum",
                        comparisonFun = "==", limitConstraints = 3))

    #'
    #' ### Case 6: Compositions of _N_ of Length _m_
    #'

    ## See Case 2. N.B. weak = TRUE has no effect
    compositionsGeneral(6, 3, repetition = TRUE)

    #'
    #' ### Case 7: Compositions of _N_ into Distinct Parts
    #'
    #' We must use `permuteGeneral` here.
    #'

    compositionsGeneral(6, 3)

    ## See Case 3
    permuteGeneral(6, 3,
                   constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 6)

    #'
    #' ### Case 8: Integer Compositions of _N_ into Parts of Varying Multiplicity
    #'

    ## compositions of 5 into 3 parts where each part can
    ## be used a maximum of 2 times.
    permuteGeneral(5, 3, freqs = rep(2, 5),
                   constraintFun = "sum",
                   comparisonFun = "==",
                   limitConstraints = 5)

    #'
    #' ## Efficiency Generating Partitions and Compositions
    #'
    #' With `compositionGeneral` we are able to take advantage of parallel computation. With `permuteGeneral`, the parallel options have no effect when generating compositions.
    #'

    ## compositions of 25
    system.time(compositionsGeneral(0:25, repetition = TRUE))

    compositionsCount(0:25, repetition=TRUE)

    ## Use multiple threads for greater efficiency. Generate
    ## over 16 million compositions in under a second!
    system.time(compositionsGeneral(0:25, repetition = TRUE, nThreads = 4))


    ## weak compositions of 12 usnig nThreads = 4
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


    ## General compositions with varying multiplicities
    system.time(comp25_gen <- permuteGeneral(25, 10, freqs = rep(4:8, 5),
                                             constraintFun = "sum",
                                             comparisonFun = "==",
                                             limitConstraints = 25))

    dim(comp25_gen)

    #'
    #' ## Safely Interrupt Execution with `cpp11::check_user_interrupt`
    #'
    #' Some of these operations can take some time, especially when you are in the exploratory phase and you don't have that much information about what type of solution you will obtain. For this reason, we have added the ability to interrupt execution. Under the hood, we call `cpp11::check_user_interrupt()` once every second to check if the user has requested for the process to be interrupted. Note that we only check for user interruptions when we cannot determine the number of results up front.
    #'
    #' This means that if we initiate a process that will take a long time or exhaust all of the available memory (e.g. we forget to put an upper limit on the number of results, relax the tolerance, etc.), we can simply hit `Ctrl + c`, or `esc` if using `RStudio`, to stop execution.
    #'

    set.seed(123)
    s = rnorm(1000)

    ## Oops!! We forgot to limit the output/put a loose tolerance
    ## There are as.numeric(comboCount(s, 20, T)) ~= 4.964324e+41
    ## This will either take a long long time, or all of your
    ## memory will be consumed!!!
    ##
    ## No problem... simply hit Ctrl + c or if in RStudio, hit esc
    ## or hit the "Stop" button

    ##
    ## system.time(testInterrupt <- partitionsGeneral(s, 20, TRUE, target = 0))
    ## Timing stopped at: 1.029 0.011 1.04
    ##

    #'
    #' ### Note about Interrupting Execution
    #'
    #' Generally, we encourage user to use iterators (See [Combinatorial Iterators in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/CombinatoricsIterators.html>)) as they offer greater flexibility. For example, with iterators it is easy to avoid resource consuming calls by only fetching a few results at a time.
    #'
    #' Here is an example of how to investigate difficult problems due to combinatorial explosion without fear of having to restart R.
    #'

    ## We use "s" defined above
    iter = partitionsIter(s, 20, TRUE, target = 0)

    ## Test one iteration to see if we need to relax the tolerance
    system.time(iter@nextIter())

    ## 8 seconds per iteration is a bit much... Let's loosen things
    ## a little by increasing the tolerance from sqrt(.Machine$double.eps)
    ## ~= 1.49e-8 to 1e-5.
    relaxedIter = partitionsIter(s, 20, TRUE, target = 0, tolerance = 1e-5)

    system.time(relaxedIter@nextIter())

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
