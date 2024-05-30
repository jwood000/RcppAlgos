reprex::reprex({
    #' The OP is looking for all integer partitions of the number 100 of maximal length of 25. The package `partitions` is equipped with a function exactly for this purpose called `restrictedparts`. E.g.:

    library(partitions)

    ## Keep the output tidy
    options(digits = 4)
    options(width = 90)

    ## all integer partitions of 10 of maximal length = 4
    restrictedparts(10, 4)

    #' Once all of the them have been generated, simply create a 5x5 matrix of each combinations (`restrictedparts` doesn't differentiate between `0 0 3` and `0 3 0`). The only problem is that there are so many possible combinations (`partitions::R(25, 100, TRUE) = 139620591`) the function throws an error when you call `restrictedparts(100, 25)`.

    test <- restrictedparts(100, 25)

    #' Since we can't generate them all via `restrictedparts`, we can generate them individually using `firstrestrictedpart` along with `nextrestrictedpart` like so:

    funPartition <- function(p, n) {
        mat <- matrix(nrow = 25, ncol = n)
        mat[, 1] <- p

        for (i in 2:n) {
            p <- nextrestrictedpart(p)
            mat[, i] <- p
        }

        mat
    }

    head(funPartition(firstrestrictedpart(100, 25), 5))

    #' The only problem here is that the iterators aren't as efficient due continuously copying.
    #'
    #' ## **Enter RcppAlgos**
    #'
    #' There is a faster way using the package `RcppAlgos` (I am the author). Similar to the `partitions` package, there is a function, `partitionsGeneral`, for generating all of partitions.

    library(RcppAlgos)
    ## Target is implicitly set to 100 below. For different targets, explicitly
    ## set the target parameter. E.g.:
    ##
    ##     partitionsGeneral(0:100, 25, TRUE, target = 200, upper = 10^5)
    ##
    ## Will generate the first 10^5 partitions of 200 using the vector 0:100

    matrixParts <- apply(
        partitionsGeneral(0:100, 25, repetition = TRUE, upper = 10^5),
        MARGIN = 1, \(x) matrix(x, ncol = 5), simplify = FALSE
    )

    all(sapply(matrixParts, sum) == 100)

    matrixParts[c(1, 90, 10^5)]

    #' ### A Better Approach: Iterators
    #'
    #' There are memory efficient iterators available as well for many topics in combinatorics including integer partitions (E.g. `partitionsIter`).
    #'
    #' Using iterators, we could create a helper function that could transform each result to our desired matrix.

    matFromIter <- function(it, ncol = 5L) {
        matrix(it@nextIter(), ncol = ncol)
    }

    ## Initialize partitions iterator
    it <- partitionsIter(0:100, 25, repetition = TRUE)
    ## Get the first 3 results
    replicate(3, matFromIter(it))
    ## Get 2 more picking up where we left off above
    replicate(2, matFromIter(it))
    ## Reset iterator
    it@startOver()
    ## Get random lexicographical result using the method: `[[`
    matrix(it[[1e6]], ncol = 5)
    ## Get the last one
    matrix(it@back(), ncol = 5)

    #' ### Need Permutations?
    #'
    #' If you really want permutations, no problem, just call `compositionsGeneral`:

    matrixComps <- apply(
        compositionsGeneral(0:100, 25, repetition = TRUE, upper = 10^5),
        MARGIN = 1, \(x) matrix(x, ncol = 5), simplify = FALSE
    )

    all(sapply(matrixComps, sum) == 100)

    ## Compare to the output of matrixCombs
    matrixComps[c(1, 90, 10^5)]

    #' ### Random Sampling
    #'
    #' Since the number of results is so massive, sampling may be our best option. Consider how many total results we are dealing with:

    partitionsCount(0:100, 25, TRUE)

    compositionsCount(0:100, 25, TRUE)

    #' We can use either `partitionsSample` or `compositionsSample` to quickly generate a candidates that can be transformed into the desired matrix output.

    ## Optional, use seed parameter for reproducibility
    apply(partitionsSample(0:100, 25, TRUE, n = 3, seed = 42), 1, \(x) {
        matrix(x, ncol = 5)
    }, simplify = FALSE)

    apply(compositionsSample(0:100, 25, TRUE, n = 3, seed = 28), 1, \(x) {
        matrix(x, ncol = 5)
    }, simplify = FALSE)

    #' ### Efficiency
    #'
    #' All of the functions are written in `C++` for ultimate efficiency. Consider iterating over 10,000 partitions.

    library(microbenchmark)
    pkg_partitions <- function(n, k, total) {
        a <- firstrestrictedpart(n, k)
        for (i in 1:(total - 1)) a <- nextrestrictedpart(a)
    }

    pkg_RcppAlgos <- function(n, k, total) {
        a <- partitionsIter(0:n, k, repetition = TRUE)
        for (i in 1:total) a@nextIter()
    }

    microbenchmark(cbRcppAlgos  = pkg_RcppAlgos(100, 25, 10^4),
                   cbPartitions = pkg_partitions(100, 25, 10^4),
                   times = 25, unit = "relative")

    #' And generating 10^5 random samples takes no time, especially when using multiple threads:

    system.time(partitionsSample(0:100, 25, TRUE, nThreads = 6,
                                 n = 1e5, seed = 42))

    system.time(compositionsSample(0:100, 25, TRUE, nThreads = 6,
                                   n = 1e5, seed = 28))

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
