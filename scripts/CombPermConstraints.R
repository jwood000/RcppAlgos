reprex::reprex({
    #'
    #' This document covers the topic of finding combinations or permutations that meet a specific set of criteria. For example, retrieving all combinations of a vector that have a product between two bounds.
    #'
    #' ### Related articles
    #'
    #' * [Integer Partitions in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/IntegerPartitions.html>)
    #' * [Integer Compositions in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/IntegerCompositions.html>)
    #' * [Attacking Problems Related to the Subset Sum Problem](<https://jwood000.github.io/RcppAlgos/articles/SubsetSum.html>)
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
    #' ## Integer Partitions & Compositions
    #'
    #' Integer partitions and compositions are natural examples of constraint-driven combinatorial enumeration. In both cases, the goal is to find results whose elements sum to a fixed target. As with combinations versus permutations, **order does not matter for partitions**, whereas **order matters for compositions**. These differences give rise to rich combinatorial structures and highly specialized algorithms for each.
    #'
    #' For this reason, they are mentioned here only briefly. For a detailed treatment of these topics, see:
    #'
    #' * [Integer Partitions in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/IntegerPartitions.html>)
    #' * [Integer Compositions in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/IntegerCompositions.html>)
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
    ## system.time(testInterrupt <- comboGeneral(
    ##     s, 20, TRUE, constraintFun = "mean", comparisonFun = "==",
    ##     limitConstraints = 0, keepResults = TRUE
    ## ))
    ## Timing stopped at: 1.993 0.008 2
    ##

    #'
    #' ### Note about Interrupting Execution
    #'
    #' Generally, we encourage user to use iterators (See [Combinatorial Iterators in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/CombinatoricsIterators.html>)) as they offer greater flexibility. For example, with iterators it is easy to avoid resource consuming calls by only fetching a few results at a time.
    #'
    #' Here is an example of how to investigate difficult problems due to combinatorial explosion without fear of having to restart R.
    #'

    ## We use "s" defined above
    iter = comboIter(s, 20, TRUE, constraintFun = "mean",
                     comparisonFun = "==", limitConstraints = 0)

    ## Test one iteration to see if we need to relax the tolerance
    system.time(iter@nextIter())

    ## That's a bit much per iteration... Let's loosen things a little by
    ## increasing the tolerance from sqrt(.Machine$double.eps) ~= 1.49e-8
    ## to 1e-5.
    relaxedIter = comboIter(s, 20, TRUE, constraintFun = "mean",
                            comparisonFun = "==", limitConstraints = 0,
                            tolerance = 1e-5)

    system.time(relaxedIter@nextIter())

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
